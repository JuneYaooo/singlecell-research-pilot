#!/usr/bin/env python3

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread


MODULE = "trajectory"
CELLTYPE_COLUMNS = ("logic_consensus_celltype", "logic_marker_rule_celltype", "celltype")


class SkipAssumption(RuntimeError):
    """Raised when trajectory assumptions are not met for this dataset."""


def append_status(out_dir: Path, module: str, status: str, detail: str = "") -> None:
    status_file = out_dir / "tables" / "module_status.tsv"
    status_file.parent.mkdir(parents=True, exist_ok=True)
    if not status_file.exists():
        status_file.write_text("module\tstatus\tdetail\n")
    clean_detail = " ".join(str(detail).replace("\t", " ").splitlines())
    with status_file.open("a") as handle:
        handle.write(f"{module}\t{status}\t{clean_detail}\n")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run a Scanpy diffusion pseudotime replacement for the Monocle trajectory module."
    )
    parser.add_argument(
        "--out",
        required=True,
        type=Path,
        help="Phase 1 output directory, for example analysis/seurat_v5_logic_run.",
    )
    parser.add_argument(
        "--target-celltype",
        default="T_cells",
        help="Preferred annotated cell type to use when it has enough cells.",
    )
    parser.add_argument(
        "--min-cells",
        default=200,
        type=int,
        help="Minimum annotated cells required for trajectory inference.",
    )
    parser.add_argument(
        "--n-hvgs",
        default=2000,
        type=int,
        help="Maximum highly variable genes to retain when Scanpy HVG selection succeeds.",
    )
    return parser.parse_args(argv)


def require_columns(frame: pd.DataFrame, columns: set[str], path: Path) -> None:
    missing = columns.difference(frame.columns)
    if missing:
        raise ValueError(f"{path} missing required columns: {', '.join(sorted(missing))}")


def choose_celltype_column(annotation: pd.DataFrame) -> str:
    for column in CELLTYPE_COLUMNS:
        if column in annotation.columns:
            return column
    raise ValueError(
        "consensus_annotation.tsv missing a usable cell type column "
        f"({', '.join(CELLTYPE_COLUMNS)})"
    )


def read_genes(path: Path) -> list[str]:
    genes = pd.read_csv(path, sep="\t", header=None)
    if genes.empty:
        raise ValueError(f"{path} contains no genes")
    return genes.iloc[:, 0].astype(str).tolist()


def read_optional_metadata(export_dir: Path) -> pd.DataFrame:
    path = export_dir / "metadata.tsv"
    if not path.exists():
        return pd.DataFrame()
    metadata = pd.read_csv(path, sep="\t")
    if "barcode" not in metadata.columns:
        return pd.DataFrame()
    metadata = metadata.copy()
    metadata["barcode"] = metadata["barcode"].astype(str)
    return metadata


def read_barcodes(export_dir: Path, metadata: pd.DataFrame, annotation: pd.DataFrame, n_cells: int) -> list[str]:
    barcode_path = export_dir / "barcodes.tsv"
    if barcode_path.exists():
        barcodes = pd.read_csv(barcode_path, sep="\t", header=None).iloc[:, 0].astype(str).tolist()
    elif not metadata.empty and len(metadata) == n_cells:
        barcodes = metadata["barcode"].astype(str).tolist()
    elif len(annotation) == n_cells:
        barcodes = annotation["barcode"].astype(str).tolist()
    else:
        raise ValueError(
            "Unable to map count matrix columns to barcodes; expected barcodes.tsv, "
            "metadata.tsv with matching rows, or annotation rows matching matrix columns"
        )

    if len(barcodes) != n_cells:
        raise ValueError(f"Barcode count mismatch: matrix has {n_cells} cells, barcodes has {len(barcodes)}")
    if len(set(barcodes)) != len(barcodes):
        raise ValueError("Barcode list contains duplicates")
    return barcodes


def load_inputs(out_dir: Path) -> tuple[sparse.csr_matrix, list[str], pd.DataFrame, pd.DataFrame, pd.DataFrame, list[str]]:
    export_dir = out_dir / "objects" / "export"
    counts_path = export_dir / "counts.mtx.gz"
    genes_path = export_dir / "genes.tsv"
    annotation_path = out_dir / "tables" / "consensus_annotation.tsv"
    umap_path = export_dir / "umap.tsv"

    for path in (counts_path, genes_path, annotation_path, umap_path):
        if not path.exists():
            raise FileNotFoundError(f"Required input missing: {path}")

    counts = mmread(counts_path)
    if not sparse.issparse(counts):
        counts = sparse.coo_matrix(counts)
    counts = counts.tocsr()

    genes = read_genes(genes_path)
    if counts.shape[0] != len(genes):
        raise ValueError(f"Gene count mismatch: matrix has {counts.shape[0]} rows, genes.tsv has {len(genes)} rows")

    annotation = pd.read_csv(annotation_path, sep="\t")
    require_columns(annotation, {"barcode", "cluster"}, annotation_path)
    annotation = annotation.copy()
    annotation["barcode"] = annotation["barcode"].astype(str)

    umap = pd.read_csv(umap_path, sep="\t")
    require_columns(umap, {"barcode", "umap_1", "umap_2"}, umap_path)
    umap = umap.copy()
    umap["barcode"] = umap["barcode"].astype(str)

    metadata = read_optional_metadata(export_dir)
    barcodes = read_barcodes(export_dir, metadata, annotation, counts.shape[1])
    return counts, genes, annotation, metadata, umap, barcodes


def select_annotation_subset(
    annotation: pd.DataFrame,
    target_celltype: str,
    min_cells: int,
) -> tuple[str, str, pd.DataFrame]:
    celltype_column = choose_celltype_column(annotation)
    clean = annotation.dropna(subset=["barcode", "cluster", celltype_column]).copy()
    clean[celltype_column] = clean[celltype_column].astype(str)
    clean["cluster"] = clean["cluster"].astype(str)

    counts = clean[celltype_column].value_counts(dropna=True)
    if target_celltype in counts.index and int(counts[target_celltype]) >= min_cells:
        selected_celltype = target_celltype
    else:
        eligible = counts[counts >= min_cells]
        if eligible.empty:
            raise SkipAssumption(f"no annotated cell type has >= {min_cells} cells")
        selected_celltype = str(eligible.sort_values(ascending=False).index[0])

    selected = clean[clean[celltype_column] == selected_celltype].copy()
    selected = selected.drop_duplicates(subset=["barcode"], keep="first")
    cluster_count = selected["cluster"].nunique(dropna=True)
    if cluster_count < 2:
        raise SkipAssumption(
            f"selected cell type {selected_celltype} has {len(selected)} cells but only {cluster_count} cluster"
        )

    return celltype_column, selected_celltype, selected


def align_counts_to_annotation(
    gene_by_cell: sparse.csr_matrix,
    genes: list[str],
    barcodes: list[str],
    selected: pd.DataFrame,
    metadata: pd.DataFrame,
) -> tuple[sparse.csr_matrix, pd.DataFrame, list[str]]:
    barcode_to_index = {barcode: i for i, barcode in enumerate(barcodes)}
    missing = sorted(set(selected["barcode"]) - set(barcode_to_index))
    if missing:
        example = ", ".join(missing[:5])
        raise ValueError(f"{len(missing)} selected annotation barcodes are absent from count matrix: {example}")

    selected = selected.copy()
    selected["_matrix_index"] = selected["barcode"].map(barcode_to_index).astype(int)
    selected = selected.sort_values("_matrix_index").reset_index(drop=True)
    matrix = gene_by_cell[:, selected["_matrix_index"].to_numpy()].T.tocsr()

    if "nCount_RNA" not in selected.columns and not metadata.empty and "nCount_RNA" in metadata.columns:
        root_metadata = metadata[["barcode", "nCount_RNA"]].drop_duplicates("barcode", keep="first")
        selected = selected.merge(root_metadata, on="barcode", how="left")

    selected = selected.drop(columns=["_matrix_index"])
    if matrix.shape[0] != len(selected):
        raise ValueError("Internal count/annotation alignment failed")
    if matrix.shape[1] != len(genes):
        raise ValueError("Internal count/gene alignment failed")
    return matrix, selected, genes


def choose_root_cell(obs: pd.DataFrame, raw_counts: sparse.csr_matrix) -> tuple[int, str]:
    cluster_sizes = obs["cluster"].astype(str).value_counts()
    smallest_size = int(cluster_sizes.min())
    smallest_clusters = sorted(cluster_sizes[cluster_sizes == smallest_size].index.astype(str).tolist())
    root_cluster = smallest_clusters[0]
    candidate_indices = np.flatnonzero(obs["cluster"].astype(str).to_numpy() == root_cluster)

    if "nCount_RNA" in obs.columns:
        n_count = pd.to_numeric(obs["nCount_RNA"], errors="coerce").to_numpy(dtype=float)
        candidate_values = n_count[candidate_indices]
        finite = np.isfinite(candidate_values)
        if finite.any():
            finite_candidates = candidate_indices[finite]
            root_idx = int(finite_candidates[np.argmin(candidate_values[finite])])
            return root_idx, f"smallest_cluster={root_cluster}; root=min_nCount_RNA"

    raw_totals = np.asarray(raw_counts.sum(axis=1)).ravel().astype(float)
    root_idx = int(candidate_indices[np.argmin(raw_totals[candidate_indices])])
    return root_idx, f"smallest_cluster={root_cluster}; root=min_raw_counts"


def maybe_select_hvgs(adata, n_hvgs: int) -> tuple[object, str]:
    import scanpy as sc

    if n_hvgs <= 0:
        return adata, "disabled"
    if adata.n_vars < 50 or adata.n_obs < 3:
        return adata, f"skipped_small_matrix genes={adata.n_vars} cells={adata.n_obs}"

    try:
        n_top_genes = min(n_hvgs, adata.n_vars)
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, flavor="seurat")
        hvg_mask = adata.var.get("highly_variable")
        if hvg_mask is None:
            return adata, "not_available"
        hvg_count = int(np.asarray(hvg_mask).sum())
        if hvg_count < 20:
            return adata, f"too_few_hvgs={hvg_count}"
        return adata[:, hvg_mask.to_numpy()].copy(), f"selected={hvg_count}"
    except Exception as exc:  # pragma: no cover - depends on Scanpy/data versions
        return adata, f"failed={type(exc).__name__}:{exc}"


def run_scanpy_trajectory(
    matrix: sparse.csr_matrix,
    genes: list[str],
    selected: pd.DataFrame,
    celltype_column: str,
    n_hvgs: int,
) -> tuple[pd.DataFrame, str]:
    import anndata as ad
    import scanpy as sc

    if matrix.shape[0] < 3 or matrix.shape[1] < 3:
        raise SkipAssumption(f"selected matrix is too small for trajectory inference: {matrix.shape}")

    obs = selected[["barcode", celltype_column, "cluster"]].copy()
    if "nCount_RNA" in selected.columns:
        obs["nCount_RNA"] = selected["nCount_RNA"]
    obs = obs.rename(columns={celltype_column: "celltype"})
    obs.index = obs["barcode"].astype(str)
    obs["celltype"] = obs["celltype"].astype(str)
    obs["cluster"] = pd.Categorical(obs["cluster"].astype(str))

    var = pd.DataFrame(index=pd.Index(genes, name="gene"))
    adata = ad.AnnData(X=matrix.copy(), obs=obs, var=var)
    adata.var_names_make_unique()

    root_idx, root_detail = choose_root_cell(selected, matrix)

    sc.pp.normalize_total(adata, target_sum=10000)
    sc.pp.log1p(adata)
    adata, hvg_detail = maybe_select_hvgs(adata, n_hvgs=n_hvgs)

    n_comps = min(50, adata.n_obs - 1, adata.n_vars - 1)
    if n_comps < 2:
        raise SkipAssumption(f"not enough dimensions for PCA after filtering: n_comps={n_comps}")
    sc.pp.pca(adata, n_comps=n_comps, svd_solver="arpack")

    n_neighbors = min(max(15, int(np.sqrt(adata.n_obs))), adata.n_obs - 1)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_comps)

    paga_detail = "not_run"
    try:
        sc.tl.paga(adata, groups="cluster")
        paga_detail = "computed"
    except Exception as exc:  # pragma: no cover - optional dependency/version behavior
        paga_detail = f"skipped={type(exc).__name__}:{exc}"

    sc.tl.diffmap(adata, n_comps=min(15, adata.n_obs - 1))
    adata.uns["iroot"] = int(root_idx)
    n_dcs = min(10, int(adata.obsm["X_diffmap"].shape[1]))
    if n_dcs < 2:
        raise SkipAssumption(f"not enough diffusion components for DPT: n_dcs={n_dcs}")
    sc.tl.dpt(adata, n_dcs=n_dcs)

    pseudotime = pd.to_numeric(adata.obs["dpt_pseudotime"], errors="coerce").to_numpy(dtype=float)
    pseudotime[~np.isfinite(pseudotime)] = np.nan
    if not np.isfinite(pseudotime).any():
        raise RuntimeError("DPT returned no finite pseudotime values")

    result = pd.DataFrame(
        {
            "barcode": adata.obs["barcode"].astype(str).to_numpy(),
            "celltype": adata.obs["celltype"].astype(str).to_numpy(),
            "cluster": adata.obs["cluster"].astype(str).to_numpy(),
            "pseudotime": pseudotime,
        }
    )
    detail = (
        f"cells={adata.n_obs}; genes={adata.n_vars}; clusters={result['cluster'].nunique()}; "
        f"neighbors={n_neighbors}; hvg={hvg_detail}; paga={paga_detail}; {root_detail}"
    )
    return result, detail


def write_pseudotime_table(result: pd.DataFrame, out_dir: Path) -> Path:
    path = out_dir / "tables" / "trajectory_pseudotime.tsv"
    result.to_csv(path, sep="\t", index=False)
    return path


def write_skip_outputs(out_dir: Path, detail: str) -> None:
    empty = pd.DataFrame(columns=["barcode", "celltype", "cluster", "pseudotime"])
    write_pseudotime_table(empty, out_dir)

    fig_path = out_dir / "figures" / "trajectory_pseudotime_umap.png"
    plt.figure(figsize=(6, 4))
    plt.axis("off")
    plt.text(0.5, 0.55, "Trajectory skipped", ha="center", va="center", fontsize=13)
    plt.text(0.5, 0.43, detail[:120], ha="center", va="center", fontsize=8, wrap=True)
    plt.tight_layout()
    plt.savefig(fig_path, dpi=180)
    plt.close()


def plot_pseudotime_umap(result: pd.DataFrame, umap: pd.DataFrame, out_dir: Path) -> Path:
    plot_df = umap.merge(result, on="barcode", how="inner")
    plot_df["umap_1"] = pd.to_numeric(plot_df["umap_1"], errors="coerce")
    plot_df["umap_2"] = pd.to_numeric(plot_df["umap_2"], errors="coerce")
    plot_df["pseudotime"] = pd.to_numeric(plot_df["pseudotime"], errors="coerce")
    plot_df = plot_df.dropna(subset=["umap_1", "umap_2", "pseudotime"])
    if plot_df.empty:
        raise RuntimeError("No finite UMAP/pseudotime rows available for trajectory plot")

    background = umap.copy()
    background["umap_1"] = pd.to_numeric(background["umap_1"], errors="coerce")
    background["umap_2"] = pd.to_numeric(background["umap_2"], errors="coerce")
    background = background.dropna(subset=["umap_1", "umap_2"])

    fig_path = out_dir / "figures" / "trajectory_pseudotime_umap.png"
    plt.figure(figsize=(6.2, 5.2))
    plt.scatter(
        background["umap_1"],
        background["umap_2"],
        c="#d9d9d9",
        s=4,
        linewidths=0,
        alpha=0.35,
    )
    points = plt.scatter(
        plot_df["umap_1"],
        plot_df["umap_2"],
        c=plot_df["pseudotime"],
        cmap="viridis",
        s=8,
        linewidths=0,
        alpha=0.9,
    )
    plt.xlabel("UMAP 1")
    plt.ylabel("UMAP 2")
    plt.title("Diffusion pseudotime")
    cbar = plt.colorbar(points)
    cbar.set_label("Pseudotime")
    plt.tight_layout()
    plt.savefig(fig_path, dpi=220)
    plt.close()
    return fig_path


def run(args: argparse.Namespace) -> str:
    out_dir = args.out
    (out_dir / "tables").mkdir(parents=True, exist_ok=True)
    (out_dir / "figures").mkdir(parents=True, exist_ok=True)

    append_status(out_dir, MODULE, "START", "scanpy_diffusion_pseudotime")
    try:
        counts, genes, annotation, metadata, umap, barcodes = load_inputs(out_dir)
        celltype_column, selected_celltype, selected = select_annotation_subset(
            annotation=annotation,
            target_celltype=args.target_celltype,
            min_cells=args.min_cells,
        )
        matrix, selected, genes = align_counts_to_annotation(
            gene_by_cell=counts,
            genes=genes,
            barcodes=barcodes,
            selected=selected,
            metadata=metadata,
        )
        result, detail = run_scanpy_trajectory(
            matrix=matrix,
            genes=genes,
            selected=selected,
            celltype_column=celltype_column,
            n_hvgs=args.n_hvgs,
        )
        write_pseudotime_table(result, out_dir)
        plot_pseudotime_umap(result, umap, out_dir)
        detail = f"celltype={selected_celltype}; {detail}"
        append_status(out_dir, MODULE, "PASS", detail)
        return detail
    except SkipAssumption as exc:
        detail = str(exc)
        write_skip_outputs(out_dir, detail)
        append_status(out_dir, MODULE, "SKIPPED_ASSUMPTION", detail)
        return detail
    except Exception as exc:
        append_status(out_dir, MODULE, "FAIL", f"{type(exc).__name__}:{exc}")
        raise


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    detail = run(args)
    print(detail)
    return 0


if __name__ == "__main__":
    sys.exit(main())
