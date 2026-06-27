#!/usr/bin/env python3

from __future__ import annotations

import argparse
import importlib.util
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread


MODULE = "cytotrace"
OFFICIAL_PACKAGE = "cytotrace2_py"
CELLTYPE_COLUMNS = (
    "logic_consensus_celltype",
    "consensus_celltype",
    "celltype",
    "annotation",
    "logic_marker_rule_celltype",
)


@dataclass(frozen=True)
class InputPaths:
    counts: Path
    genes: Path
    barcodes: Path
    annotations: Path
    umap: Path


def append_status(out_dir: Path, module: str, status: str, detail: str = "") -> None:
    status_file = out_dir / "tables" / "module_status.tsv"
    status_file.parent.mkdir(parents=True, exist_ok=True)
    if not status_file.exists():
        status_file.write_text("module\tstatus\tdetail\n")
    clean_detail = " ".join(str(detail).replace("\t", " ").splitlines())
    with status_file.open("a") as handle:
        handle.write(f"{module}\t{status}\t{clean_detail}\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run a bounded CytoTRACE-like gene-count potential module. "
            "This produces the user-facing CytoTRACE workflow outputs, but "
            "does not claim equivalence to the official CytoTRACE/CytoTRACE2 package."
        )
    )
    parser.add_argument("--out", required=True, type=Path, help="Logic output directory.")
    parser.add_argument("--counts", type=Path, help="MatrixMarket counts path. Defaults under --out.")
    parser.add_argument("--genes", type=Path, help="genes.tsv path. Defaults under --out.")
    parser.add_argument("--barcodes", type=Path, help="barcodes.tsv path. Defaults under --out.")
    parser.add_argument("--annotations", type=Path, help="consensus_annotation.tsv path. Defaults under --out.")
    parser.add_argument("--umap", type=Path, help="UMAP TSV path. Defaults under --out.")
    parser.add_argument(
        "--run-official-cytotrace2",
        action="store_true",
        help=(
            "Attempt the official cytotrace2-py runner. The default is detection-only "
            "because the official API can require dense input materialization and model download."
        ),
    )
    parser.add_argument(
        "--official-species",
        default="human",
        choices=("human", "mouse"),
        help="Species argument to use if --run-official-cytotrace2 is enabled.",
    )
    return parser.parse_args()


def default_paths(out_dir: Path) -> InputPaths:
    export_dir = out_dir / "objects" / "export"
    consensus_umap = out_dir / "tables" / "annotated_consensus_umap.tsv"
    return InputPaths(
        counts=export_dir / "counts.mtx.gz",
        genes=export_dir / "genes.tsv",
        barcodes=export_dir / "barcodes.tsv",
        annotations=out_dir / "tables" / "consensus_annotation.tsv",
        umap=consensus_umap if consensus_umap.exists() else export_dir / "umap.tsv",
    )


def resolve_paths(args: argparse.Namespace) -> InputPaths:
    defaults = default_paths(args.out)
    return InputPaths(
        counts=args.counts or defaults.counts,
        genes=args.genes or defaults.genes,
        barcodes=args.barcodes or defaults.barcodes,
        annotations=args.annotations or defaults.annotations,
        umap=args.umap or defaults.umap,
    )


def require_files(paths: InputPaths) -> None:
    missing = [
        str(path)
        for path in (paths.counts, paths.genes, paths.barcodes, paths.annotations, paths.umap)
        if not path.exists()
    ]
    if missing:
        raise FileNotFoundError("Missing required input file(s): " + ", ".join(missing))


def read_vector(path: Path, name: str) -> list[str]:
    table = pd.read_csv(path, sep="\t", header=None, dtype=str)
    if table.empty:
        raise ValueError(f"{name} is empty: {path}")
    values = table.iloc[:, 0].fillna("").astype(str).str.strip().tolist()
    if not values:
        raise ValueError(f"{name} has no values: {path}")
    return values


def choose_celltype_column(annotations: pd.DataFrame) -> str:
    for column in CELLTYPE_COLUMNS:
        if column in annotations.columns:
            return column
    raise ValueError(
        "consensus_annotation.tsv missing a usable cell type column "
        f"({', '.join(CELLTYPE_COLUMNS)})"
    )


def read_annotations(path: Path) -> tuple[pd.DataFrame, str]:
    annotations = pd.read_csv(path, sep="\t", dtype=str)
    if annotations.empty:
        raise ValueError(f"Annotation file is empty: {path}")
    if "barcode" not in annotations.columns:
        raise ValueError(f"{path} must include a barcode column")

    celltype_column = choose_celltype_column(annotations)
    keep = ["barcode", celltype_column]
    keep.extend(column for column in ("sample_id", "cluster") if column in annotations.columns)
    clean = annotations[keep].rename(columns={celltype_column: "celltype"}).copy()
    clean["barcode"] = clean["barcode"].astype(str)
    clean["celltype"] = clean["celltype"].astype(str)
    clean = clean.drop_duplicates(subset=["barcode"], keep="first")
    return clean, celltype_column


def read_counts(path: Path, genes: list[str], barcodes: list[str]) -> sparse.csr_matrix:
    matrix = mmread(str(path))
    if not sparse.issparse(matrix):
        matrix = sparse.coo_matrix(matrix)
    matrix = matrix.tocsr()

    if matrix.shape == (len(genes), len(barcodes)):
        return matrix
    if matrix.shape == (len(barcodes), len(genes)):
        return matrix.T.tocsr()
    raise ValueError(
        "Counts matrix shape does not match genes/barcodes: "
        f"matrix={matrix.shape}, genes={len(genes)}, barcodes={len(barcodes)}"
    )


def align_to_annotations(
    gene_by_cell: sparse.csr_matrix,
    annotations: pd.DataFrame,
    barcodes: list[str],
) -> tuple[sparse.csr_matrix, pd.DataFrame]:
    barcode_index = pd.Series(np.arange(len(barcodes)), index=barcodes)
    aligned = annotations.copy()
    aligned["_matrix_col"] = aligned["barcode"].map(barcode_index)
    aligned = aligned.dropna(subset=["_matrix_col"]).copy()
    if aligned.empty:
        raise ValueError("No annotated barcodes were found in the exported count matrix")
    aligned["_matrix_col"] = aligned["_matrix_col"].astype(int)
    aligned = aligned.sort_values("_matrix_col").reset_index(drop=True)
    matrix = gene_by_cell[:, aligned["_matrix_col"].to_numpy()].tocsr()
    aligned = aligned.drop(columns=["_matrix_col"])
    return matrix, aligned


def percentile_rank(values: np.ndarray) -> np.ndarray:
    series = pd.Series(values)
    if len(series) <= 1:
        return np.ones(len(series), dtype=float)
    ranks = series.rank(method="average", pct=True).to_numpy(dtype=float)
    return np.clip(ranks, 0.0, 1.0)


def score_cytotrace_like(gene_by_cell: sparse.csr_matrix, metadata: pd.DataFrame) -> pd.DataFrame:
    detected = np.asarray(gene_by_cell.getnnz(axis=0)).ravel().astype(float)
    total_counts = np.asarray(gene_by_cell.sum(axis=0)).ravel().astype(float)
    detected_rank = percentile_rank(detected)
    total_rank = percentile_rank(np.log1p(total_counts))

    scores = metadata.copy()
    scores["n_detected_genes"] = detected.astype(int)
    scores["total_counts"] = total_counts
    scores["cytotrace_like_score"] = (0.8 * detected_rank) + (0.2 * total_rank)
    scores["cytotrace_like_score"] = scores["cytotrace_like_score"].clip(0.0, 1.0)
    scores["method"] = "gene_count_rank_0.8_detected_0.2_total_counts"
    scores["equivalence_level"] = "TRACEABLE_REPLACEMENT_NOT_OFFICIAL_CYTOTRACE"
    return scores


def method_status_row(
    method: str,
    status: str,
    detail: str,
    output_table: str,
    *,
    official_package_available: bool | None = None,
    official_execution: str = "",
) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "module": MODULE,
                "method": method,
                "status": status,
                "detail": detail,
                "output_table": output_table,
                "official_package": OFFICIAL_PACKAGE,
                "official_package_available": "" if official_package_available is None else str(official_package_available),
                "official_execution": official_execution,
            }
        ]
    )


def try_official_cytotrace2(
    gene_by_cell: sparse.csr_matrix,
    metadata: pd.DataFrame,
    genes: list[str],
    *,
    run_official: bool,
    species: str,
) -> tuple[pd.DataFrame | None, str]:
    """Try a bounded official CytoTRACE2 Python execution when the package API is available.

    CytoTRACE2 package APIs have changed across distributions. This function only runs
    when a recognizable callable is present; otherwise the workflow falls back to the
    deterministic local score and records the reason.
    """
    if not package_available(OFFICIAL_PACKAGE):
        return None, f"{OFFICIAL_PACKAGE} package not installed"

    if not run_official:
        return (
            None,
            (
                f"{OFFICIAL_PACKAGE} package installed; official execution disabled by default "
                "because it requires dense file input and may download model weights; rerun "
                "with --run-official-cytotrace2 to attempt official CytoTRACE2"
            ),
        )

    try:
        from cytotrace2_py import cytotrace2_py as cytotrace2  # type: ignore
    except Exception as exc:  # pragma: no cover - depends on optional dependency
        return None, f"{OFFICIAL_PACKAGE} import failed: {type(exc).__name__}:{exc}"

    callable_names = ("cytotrace2", "predict", "run_cytotrace2", "CytoTRACE2")
    runner = next((getattr(cytotrace2, name, None) for name in callable_names if callable(getattr(cytotrace2, name, None))), None)
    if runner is None:
        exposed = ",".join(sorted(name for name in dir(cytotrace2) if not name.startswith("_"))[:20])
        return None, f"cytotrace2 import succeeded but no supported callable found; exposed={exposed}"

    try:  # pragma: no cover - optional package not present in validated env
        cell_by_gene = gene_by_cell.T.tocsr().astype(np.float32)
        expression = pd.DataFrame.sparse.from_spmatrix(cell_by_gene, index=metadata["barcode"], columns=genes)
        try:
            result = runner(expression, species=species)
        except TypeError:
            result = runner(expression)
        if isinstance(result, pd.DataFrame):
            table = result.copy()
        else:
            table = pd.DataFrame(result)
        if "barcode" not in table.columns:
            table.insert(0, "barcode", metadata["barcode"].to_numpy())
        table = metadata.merge(table, on="barcode", how="left")
        table["n_detected_genes"] = np.asarray(gene_by_cell.getnnz(axis=0)).ravel().astype(int)
        table["total_counts"] = np.asarray(gene_by_cell.sum(axis=0)).ravel().astype(float)
        detected_rank = percentile_rank(table["n_detected_genes"].to_numpy(dtype=float))
        total_rank = percentile_rank(np.log1p(table["total_counts"].to_numpy(dtype=float)))
        table["cytotrace_like_score"] = np.clip((0.8 * detected_rank) + (0.2 * total_rank), 0.0, 1.0)
        if "cytotrace2_score" not in table.columns:
            numeric_candidates = [
                column
                for column in table.columns
                if column not in {"barcode", "celltype", "sample_id", "cluster"}
                and pd.api.types.is_numeric_dtype(table[column])
            ]
            if numeric_candidates:
                table["cytotrace2_score"] = pd.to_numeric(table[numeric_candidates[0]], errors="coerce")
        table["method"] = "official_cytotrace2"
        table["equivalence_level"] = "OFFICIAL_CYTOTRACE2"
        return table, "official CytoTRACE2 callable executed"
    except Exception as exc:
        return None, f"cytotrace2 execution failed: {type(exc).__name__}:{exc}"


def summarize_scores(scores: pd.DataFrame) -> pd.DataFrame:
    score_column = "cytotrace2_score" if "cytotrace2_score" in scores.columns else "cytotrace_like_score"
    grouped = scores.groupby("celltype", dropna=False)
    summary = grouped.agg(
        cells=("barcode", "count"),
        mean_score=(score_column, "mean"),
        median_score=(score_column, "median"),
        mean_detected_genes=("n_detected_genes", "mean"),
        median_detected_genes=("n_detected_genes", "median"),
    ).reset_index()
    return summary.sort_values(["median_score", "cells"], ascending=[False, False])


def plot_umap(scores: pd.DataFrame, umap_path: Path, out_dir: Path) -> None:
    score_column = "cytotrace2_score" if "cytotrace2_score" in scores.columns else "cytotrace_like_score"
    umap = pd.read_csv(umap_path, sep="\t", dtype=str)
    required = {"barcode", "umap_1", "umap_2"}
    missing = required.difference(umap.columns)
    if missing:
        raise ValueError(f"{umap_path} missing required columns: {', '.join(sorted(missing))}")

    plot_df = umap[["barcode", "umap_1", "umap_2"]].merge(
        scores[["barcode", score_column]], on="barcode", how="left"
    )
    for column in ("umap_1", "umap_2", score_column):
        plot_df[column] = pd.to_numeric(plot_df[column], errors="coerce")
    plot_df = plot_df.dropna(subset=["umap_1", "umap_2"])
    scored = plot_df.dropna(subset=[score_column])

    plt.figure(figsize=(6.5, 5.5))
    plt.scatter(plot_df["umap_1"], plot_df["umap_2"], c="#d9d9d9", s=3, linewidths=0, alpha=0.45)
    if not scored.empty:
        points = plt.scatter(
            scored["umap_1"],
            scored["umap_2"],
            c=scored[score_column],
            cmap="viridis",
            s=5,
            linewidths=0,
            alpha=0.9,
            vmin=0,
            vmax=1,
        )
        cbar = plt.colorbar(points, fraction=0.046, pad=0.04)
        cbar.set_label("CytoTRACE-like potential score")
    plt.xlabel("UMAP 1")
    plt.ylabel("UMAP 2")
    plt.title("CytoTRACE-like differentiation potential")
    plt.tight_layout()
    plt.savefig(out_dir / "figures" / "cytotrace_umap.png", dpi=220)
    plt.close()


def package_available(name: str) -> bool:
    return importlib.util.find_spec(name) is not None


def run(args: argparse.Namespace) -> str:
    out_dir = args.out
    (out_dir / "tables").mkdir(parents=True, exist_ok=True)
    (out_dir / "figures").mkdir(parents=True, exist_ok=True)
    append_status(out_dir, MODULE, "START", "cytotrace_like_gene_count_potential")

    paths = resolve_paths(args)
    require_files(paths)
    genes = read_vector(paths.genes, "genes.tsv")
    barcodes = read_vector(paths.barcodes, "barcodes.tsv")
    annotations, celltype_column = read_annotations(paths.annotations)
    gene_by_cell = read_counts(paths.counts, genes, barcodes)
    aligned_counts, aligned_annotations = align_to_annotations(gene_by_cell, annotations, barcodes)

    official_available = package_available(OFFICIAL_PACKAGE)
    official_scores, official_detail = try_official_cytotrace2(
        aligned_counts,
        aligned_annotations,
        genes,
        run_official=args.run_official_cytotrace2,
        species=args.official_species,
    )
    if official_scores is not None:
        scores = official_scores
        method_status = method_status_row(
            "official_cytotrace2",
            "PASS",
            official_detail,
            "tables/cytotrace_scores.tsv",
            official_package_available=official_available,
            official_execution="attempted",
        )
    else:
        scores = score_cytotrace_like(aligned_counts, aligned_annotations)
        if not official_available:
            official_execution = "not_available"
        elif not args.run_official_cytotrace2:
            official_execution = "disabled_by_default"
        else:
            official_execution = "attempted_failed"
        method_status = method_status_row(
            "cytotrace_like_gene_count_rank",
            "FALLBACK",
            official_detail,
            "tables/cytotrace_scores.tsv",
            official_package_available=official_available,
            official_execution=official_execution,
        )
    summary = summarize_scores(scores)
    scores.to_csv(out_dir / "tables" / "cytotrace_scores.tsv", sep="\t", index=False)
    summary.to_csv(out_dir / "tables" / "cytotrace_celltype_summary.tsv", sep="\t", index=False)
    method_status.to_csv(out_dir / "tables" / "cytotrace_method_status.tsv", sep="\t", index=False)
    plot_umap(scores, paths.umap, out_dir)

    detail = (
        f"cells={len(scores)} celltypes={scores['celltype'].nunique()} "
        f"celltype_column={celltype_column} method={method_status.iloc[0]['method']} "
        f"method_status={method_status.iloc[0]['status']} official_detail={official_detail}"
    )
    append_status(out_dir, MODULE, "PASS", detail)
    return detail


def main() -> None:
    args = parse_args()
    try:
        detail = run(args)
    except Exception as exc:
        try:
            append_status(args.out, MODULE, "FAIL", f"{type(exc).__name__}:{exc}")
        except Exception:
            pass
        raise
    print(detail)


if __name__ == "__main__":
    main()
