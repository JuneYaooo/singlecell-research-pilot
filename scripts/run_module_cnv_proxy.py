#!/usr/bin/env python3

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


MODULE = "cnv_proxy"
TUMOR_CELLTYPE = "Epithelial_tumor_like"
REFERENCE_CELLTYPES = {
    "T_cells",
    "B_cells",
    "Plasma",
    "Myeloid_monocyte",
    "Endothelial",
    "Fibroblast_stromal",
    "Mast",
}
CELLTYPE_COLUMNS = (
    "logic_consensus_celltype",
    "celltype",
    "annotation",
    "logic_marker_rule_celltype",
)


@dataclass(frozen=True)
class InputPaths:
    counts: Path
    genes: Path
    annotations: Path
    umap: Path
    barcodes: Path


def append_status(out_dir: Path, module: str, status: str, detail: str = "") -> None:
    status_file = out_dir / "tables" / "module_status.tsv"
    status_file.parent.mkdir(parents=True, exist_ok=True)
    if not status_file.exists():
        status_file.write_text("module\tstatus\tdetail\n")
    clean_detail = " ".join(str(detail).replace("\t", " ").splitlines())
    with status_file.open("a") as handle:
        handle.write(f"{module}\t{status}\t{clean_detail}\n")


def default_paths(out_dir: Path) -> InputPaths:
    export_dir = out_dir / "objects" / "export"
    return InputPaths(
        counts=export_dir / "counts.mtx.gz",
        genes=export_dir / "genes.tsv",
        annotations=out_dir / "tables" / "consensus_annotation.tsv",
        umap=export_dir / "umap.tsv",
        barcodes=export_dir / "barcodes.tsv",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compute a bounded expression-deviation CNV proxy. This is not "
            "copykat, inferCNV, or a chromosome-ordered CNV caller."
        )
    )
    parser.add_argument("--out", required=True, type=Path, help="core workflow output directory.")
    parser.add_argument("--counts", type=Path, help="MatrixMarket counts path. Defaults under --out.")
    parser.add_argument("--genes", type=Path, help="Gene name TSV path. Defaults under --out.")
    parser.add_argument("--annotations", type=Path, help="Consensus annotation TSV path. Defaults under --out.")
    parser.add_argument("--umap", type=Path, help="UMAP TSV path. Defaults under --out.")
    parser.add_argument("--barcodes", type=Path, help="Barcode TSV path used to align matrix columns. Defaults under --out.")
    parser.add_argument("--top-genes", type=int, default=2000, help="Number of variable genes to score.")
    parser.add_argument("--normalization-scale", type=float, default=10000.0, help="Per-cell count normalization scale.")
    parser.add_argument(
        "--min-reference-cells",
        type=int,
        default=10,
        help="Minimum reference cells required to estimate mean and standard deviation.",
    )
    parser.add_argument(
        "--min-epithelial-cells",
        type=int,
        default=1,
        help="Minimum Epithelial_tumor_like cells required to score.",
    )
    return parser.parse_args()


def resolve_paths(args: argparse.Namespace) -> InputPaths:
    defaults = default_paths(args.out)
    return InputPaths(
        counts=args.counts or defaults.counts,
        genes=args.genes or defaults.genes,
        annotations=args.annotations or defaults.annotations,
        umap=args.umap or defaults.umap,
        barcodes=args.barcodes or defaults.barcodes,
    )


def require_files(paths: InputPaths) -> None:
    missing = [
        str(path)
        for path in (paths.counts, paths.genes, paths.annotations, paths.umap, paths.barcodes)
        if not path.exists()
    ]
    if missing:
        raise FileNotFoundError("Missing required input file(s): " + ", ".join(missing))


def read_vector(path: Path, name: str) -> list[str]:
    values = pd.read_csv(path, sep="\t", header=None, dtype=str).iloc[:, 0].fillna("").tolist()
    if not values:
        raise ValueError(f"{name} is empty: {path}")
    return values


def read_annotations(path: Path) -> pd.DataFrame:
    annotations = pd.read_csv(path, sep="\t", dtype=str)
    if "barcode" not in annotations.columns:
        raise ValueError(f"{path} must include a barcode column")

    celltype_column = next((col for col in CELLTYPE_COLUMNS if col in annotations.columns), None)
    if celltype_column is None:
        raise ValueError(
            f"{path} must include one of these cell type columns: {', '.join(CELLTYPE_COLUMNS)}"
        )

    clean = annotations[["barcode", celltype_column]].rename(columns={celltype_column: "celltype"}).copy()
    clean["barcode"] = clean["barcode"].astype(str)
    clean["celltype"] = clean["celltype"].astype(str)
    clean = clean.drop_duplicates(subset=["barcode"], keep="first")
    return clean


def load_counts_matrix(path: Path, genes: list[str], barcodes: list[str]) -> sparse.csr_matrix:
    matrix = mmread(path).tocsr()
    if matrix.shape == (len(genes), len(barcodes)):
        matrix = matrix.T.tocsr()
    elif matrix.shape == (len(barcodes), len(genes)):
        matrix = matrix.tocsr()
    else:
        raise ValueError(
            "Counts matrix shape does not match genes/barcodes: "
            f"matrix={matrix.shape}, genes={len(genes)}, barcodes={len(barcodes)}"
        )
    return matrix


def normalize_log1p(matrix: sparse.csr_matrix, scale: float) -> sparse.csr_matrix:
    if scale <= 0:
        raise ValueError("--normalization-scale must be positive")

    totals = np.asarray(matrix.sum(axis=1)).ravel().astype(float)
    factors = np.zeros_like(totals, dtype=float)
    nonzero = totals > 0
    factors[nonzero] = scale / totals[nonzero]

    normalized = matrix.multiply(factors[:, None]).tocsr()
    normalized.data = np.log1p(normalized.data)
    return normalized


def select_variable_genes(matrix: sparse.csr_matrix, genes: list[str], top_n: int) -> np.ndarray:
    if top_n <= 0:
        raise ValueError("--top-genes must be positive")

    means = np.asarray(matrix.mean(axis=0)).ravel()
    squared_means = np.asarray(matrix.power(2).mean(axis=0)).ravel()
    variances = squared_means - np.square(means)
    valid = np.flatnonzero(np.isfinite(variances) & (variances > 0))
    if valid.size == 0:
        raise RuntimeError("No variable genes available after normalization")

    n_select = min(int(top_n), valid.size, len(genes))
    order = valid[np.argsort(variances[valid], kind="mergesort")[-n_select:]]
    return order[np.argsort(variances[order])[::-1]]


def score_expression_deviation(
    normalized: sparse.csr_matrix,
    metadata: pd.DataFrame,
    genes: list[str],
    top_genes: int,
) -> tuple[pd.DataFrame, int]:
    gene_idx = select_variable_genes(normalized, genes, top_genes)
    dense = normalized[:, gene_idx].toarray().astype(np.float32, copy=False)

    reference_mask = metadata["celltype"].isin(REFERENCE_CELLTYPES).to_numpy()
    reference_values = dense[reference_mask, :]
    reference_mean = reference_values.mean(axis=0)
    reference_sd = reference_values.std(axis=0)

    usable = np.isfinite(reference_mean) & np.isfinite(reference_sd) & (reference_sd > 0)
    if not np.any(usable):
        raise RuntimeError("No selected genes have nonzero reference standard deviation")

    z_scores = (dense[:, usable] - reference_mean[usable]) / reference_sd[usable]
    result = metadata[["barcode", "celltype"]].copy()
    result["cnv_proxy_score"] = np.mean(np.abs(z_scores), axis=1)
    result["cnv_proxy_signed_score"] = np.mean(z_scores, axis=1)
    return result, int(np.sum(usable))


def assumption_failure(metadata: pd.DataFrame, min_epithelial: int, min_reference: int) -> str | None:
    tumor_count = int((metadata["celltype"] == TUMOR_CELLTYPE).sum())
    reference_count = int(metadata["celltype"].isin(REFERENCE_CELLTYPES).sum())
    if tumor_count < min_epithelial:
        return f"requires {TUMOR_CELLTYPE} cells; found={tumor_count}"
    if reference_count < min_reference:
        return (
            "requires reference cells from "
            f"{','.join(sorted(REFERENCE_CELLTYPES))}; found={reference_count}"
        )
    return None


def plot_umap(scores: pd.DataFrame, umap: pd.DataFrame, out_dir: Path) -> None:
    required = {"barcode", "umap_1", "umap_2"}
    missing = required.difference(umap.columns)
    if missing:
        raise ValueError(f"UMAP table missing required columns: {', '.join(sorted(missing))}")

    plot_df = umap[["barcode", "umap_1", "umap_2"]].merge(scores, on="barcode", how="left")
    plot_df["umap_1"] = pd.to_numeric(plot_df["umap_1"], errors="coerce")
    plot_df["umap_2"] = pd.to_numeric(plot_df["umap_2"], errors="coerce")
    plot_df["cnv_proxy_score"] = pd.to_numeric(plot_df["cnv_proxy_score"], errors="coerce")
    plot_df = plot_df.dropna(subset=["umap_1", "umap_2"])

    scored = plot_df.dropna(subset=["cnv_proxy_score"])
    plt.figure(figsize=(6.5, 5.5))
    plt.scatter(plot_df["umap_1"], plot_df["umap_2"], c="#d9d9d9", s=3, linewidths=0, alpha=0.45)
    if not scored.empty:
        points = plt.scatter(
            scored["umap_1"],
            scored["umap_2"],
            c=scored["cnv_proxy_score"],
            cmap="magma",
            s=5,
            linewidths=0,
            alpha=0.85,
        )
        cbar = plt.colorbar(points, fraction=0.046, pad=0.04)
        cbar.set_label("CNV proxy score")
    plt.xlabel("UMAP 1")
    plt.ylabel("UMAP 2")
    plt.title("Expression-deviation CNV proxy")
    plt.tight_layout()
    plt.savefig(out_dir / "figures" / "cnv_proxy_umap.png", dpi=220)
    plt.close()


def run(args: argparse.Namespace) -> str:
    out_dir = args.out
    (out_dir / "tables").mkdir(parents=True, exist_ok=True)
    (out_dir / "figures").mkdir(parents=True, exist_ok=True)

    paths = resolve_paths(args)
    append_status(out_dir, MODULE, "START", "expression_deviation_cnv_proxy")

    require_files(paths)
    genes = read_vector(paths.genes, "genes.tsv")
    barcodes = read_vector(paths.barcodes, "barcodes.tsv")
    annotations = read_annotations(paths.annotations)
    matrix = load_counts_matrix(paths.counts, genes, barcodes)

    matrix_metadata = pd.DataFrame({"barcode": barcodes})
    metadata = matrix_metadata.merge(annotations, on="barcode", how="inner")
    metadata = metadata[metadata["celltype"].isin({TUMOR_CELLTYPE, *REFERENCE_CELLTYPES})].copy()

    assumption_detail = assumption_failure(metadata, args.min_epithelial_cells, args.min_reference_cells)
    if assumption_detail is not None:
        append_status(out_dir, MODULE, "SKIPPED_ASSUMPTION", assumption_detail)
        return assumption_detail

    row_lookup = pd.Series(np.arange(len(barcodes)), index=barcodes)
    row_indices = row_lookup.loc[metadata["barcode"]].to_numpy()
    normalized = normalize_log1p(matrix[row_indices, :], args.normalization_scale)
    scores, usable_genes = score_expression_deviation(normalized, metadata, genes, args.top_genes)

    scores.to_csv(out_dir / "tables" / "cnv_proxy_scores.tsv", sep="\t", index=False)
    umap = pd.read_csv(paths.umap, sep="\t", dtype=str)
    plot_umap(scores, umap, out_dir)

    tumor_count = int((metadata["celltype"] == TUMOR_CELLTYPE).sum())
    reference_count = int(metadata["celltype"].isin(REFERENCE_CELLTYPES).sum())
    detail = (
        f"scored_cells={len(scores)} epithelial_cells={tumor_count} "
        f"reference_cells={reference_count} variable_genes={usable_genes}"
    )
    append_status(out_dir, MODULE, "PASS", detail)
    return detail


def main() -> None:
    args = parse_args()
    try:
        detail = run(args)
    except Exception as exc:
        out_dir = args.out
        try:
            append_status(out_dir, MODULE, "FAIL", f"{type(exc).__name__}:{exc}")
        except Exception:
            pass
        raise
    print(detail)


if __name__ == "__main__":
    main()
