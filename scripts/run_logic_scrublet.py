#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy import sparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def append_status(out_dir: Path, module: str, status: str, detail: str = "") -> None:
    status_file = out_dir / "tables" / "module_status.tsv"
    status_file.parent.mkdir(parents=True, exist_ok=True)
    if not status_file.exists():
        status_file.write_text("module\tstatus\tdetail\n")
    clean_detail = " ".join(str(detail).replace("\t", " ").splitlines())
    with status_file.open("a") as handle:
        handle.write(f"{module}\t{status}\t{clean_detail}\n")


def load_inputs(export_dir: Path) -> tuple[sparse.csr_matrix, pd.DataFrame, pd.DataFrame]:
    matrix = mmread(export_dir / "counts.mtx.gz").tocsr().T
    metadata = pd.read_csv(export_dir / "metadata.tsv", sep="\t")
    umap_path = export_dir / "umap.tsv"
    umap = pd.read_csv(umap_path, sep="\t") if umap_path.exists() else pd.DataFrame()
    if "barcode" not in metadata.columns:
        raise ValueError("metadata.tsv must include a barcode column")
    if matrix.shape[0] != metadata.shape[0]:
        raise ValueError(
            f"Cell count mismatch: matrix has {matrix.shape[0]} cells, metadata has {metadata.shape[0]} rows"
        )
    return matrix, metadata, umap


def fallback_scores(sample_matrix: sparse.csr_matrix, sample_meta: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    if {"nCount_RNA", "nFeature_RNA"}.issubset(sample_meta.columns):
        counts = sample_meta["nCount_RNA"].to_numpy(dtype=float)
        features = sample_meta["nFeature_RNA"].to_numpy(dtype=float)
    else:
        counts = np.asarray(sample_matrix.sum(axis=1)).ravel().astype(float)
        features = np.asarray((sample_matrix > 0).sum(axis=1)).ravel().astype(float)

    def robust_z(values: np.ndarray) -> np.ndarray:
        med = np.nanmedian(values)
        mad = np.nanmedian(np.abs(values - med))
        scale = mad if mad > 0 else np.nanstd(values)
        if not np.isfinite(scale) or scale == 0:
            return np.zeros_like(values, dtype=float)
        return (values - med) / scale

    scores = robust_z(counts) + robust_z(features)
    cutoff = np.nanquantile(scores, 0.95)
    predicted = scores >= cutoff
    score_min = np.nanmin(scores)
    score_max = np.nanmax(scores)
    if score_max > score_min:
        scores = (scores - score_min) / (score_max - score_min)
    else:
        scores = np.zeros_like(scores, dtype=float)
    return scores, predicted


def run_scrublet_by_sample(matrix: sparse.csr_matrix, metadata: pd.DataFrame) -> pd.DataFrame:
    try:
        import scrublet as scr

        scrublet_available = True
        scrublet_error = ""
    except Exception as exc:  # pragma: no cover - environment-dependent
        scr = None
        scrublet_available = False
        scrublet_error = repr(exc)

    sample_column = "sample_id" if "sample_id" in metadata.columns else None
    sample_values = metadata[sample_column].astype(str) if sample_column else pd.Series(["all"] * len(metadata))

    outputs: list[pd.DataFrame] = []
    for sample_id in sample_values.unique():
        idx = np.flatnonzero(sample_values.to_numpy() == sample_id)
        sample_matrix = matrix[idx, :]
        sample_meta = metadata.iloc[idx].copy()

        method = "scrublet"
        detail = ""
        if scrublet_available:
            try:
                scrub = scr.Scrublet(sample_matrix, expected_doublet_rate=0.06)
                n_prin_comps = min(30, max(2, sample_matrix.shape[0] - 1), max(2, sample_matrix.shape[1] - 1))
                scores, predicted = scrub.scrub_doublets(
                    min_counts=2,
                    min_cells=3,
                    min_gene_variability_pctl=85,
                    n_prin_comps=n_prin_comps,
                    verbose=False,
                )
                score_min = np.nanmin(scores)
                score_max = np.nanmax(scores)
                if predicted is None or not np.isfinite(score_min) or not np.isfinite(score_max) or score_max == score_min:
                    scores, predicted = fallback_scores(sample_matrix, sample_meta)
                    method = "fallback_qc_top5pct"
                    detail = "scrublet_no_threshold_or_constant_scores"
            except Exception as exc:  # pragma: no cover - data-dependent fallback
                scores, predicted = fallback_scores(sample_matrix, sample_meta)
                method = "fallback_qc_top5pct"
                detail = f"scrublet_failed:{type(exc).__name__}:{exc}"
        else:
            scores, predicted = fallback_scores(sample_matrix, sample_meta)
            method = "fallback_qc_top5pct"
            detail = f"scrublet_unavailable:{scrublet_error}"

        outputs.append(
            pd.DataFrame(
                {
                    "barcode": sample_meta["barcode"].to_numpy(),
                    "sample_id": sample_id,
                    "doublet_score": scores,
                    "predicted_doublet": predicted.astype(bool),
                    "method": method,
                    "detail": detail,
                }
            )
        )

    result = pd.concat(outputs, axis=0, ignore_index=True)
    result["predicted_doublet"] = result["predicted_doublet"].map({True: "TRUE", False: "FALSE"})
    return result


def write_counts(result: pd.DataFrame, out_dir: Path) -> None:
    counts = (
        result.assign(predicted_doublet=result["predicted_doublet"] == "TRUE")
        .groupby(["sample_id", "predicted_doublet", "method"], dropna=False)
        .size()
        .reset_index(name="cells")
    )
    before = result.groupby("sample_id", dropna=False).size().reset_index(name="cells_before")
    after = (
        result[result["predicted_doublet"] != "TRUE"]
        .groupby("sample_id", dropna=False)
        .size()
        .reset_index(name="cells_after")
    )
    summary = before.merge(after, on="sample_id", how="left")
    summary["cells_after"] = summary["cells_after"].fillna(0).astype(int)
    summary["predicted_doublets"] = summary["cells_before"] - summary["cells_after"]
    summary.to_csv(out_dir / "tables" / "cell_counts_before_after_doublet.tsv", sep="\t", index=False)
    counts.to_csv(out_dir / "tables" / "doublet_method_counts.tsv", sep="\t", index=False)


def plot_outputs(result: pd.DataFrame, umap: pd.DataFrame, out_dir: Path) -> None:
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(7, 4))
    for method, sub in result.groupby("method"):
        plt.hist(sub["doublet_score"], bins=40, alpha=0.6, label=method)
    plt.xlabel("Doublet score")
    plt.ylabel("Cells")
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(fig_dir / "doublet_score_histogram.png", dpi=200)
    plt.close()

    if {"barcode", "umap_1", "umap_2"}.issubset(umap.columns):
        plot_df = umap.merge(result[["barcode", "predicted_doublet"]], on="barcode", how="left")
        colors = plot_df["predicted_doublet"].map({"TRUE": "#d73027", "FALSE": "#4575b4"}).fillna("#999999")
        plt.figure(figsize=(6, 5))
        plt.scatter(plot_df["umap_1"], plot_df["umap_2"], c=colors, s=4, linewidths=0, alpha=0.8)
        plt.xlabel("UMAP 1")
        plt.ylabel("UMAP 2")
        plt.title("Predicted doublets")
        plt.tight_layout()
        plt.savefig(fig_dir / "umap_doublet_class.png", dpi=220)
        plt.close()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--export-dir", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    export_dir = Path(args.export_dir)
    out_dir = Path(args.out)
    (out_dir / "tables").mkdir(parents=True, exist_ok=True)
    (out_dir / "figures").mkdir(parents=True, exist_ok=True)

    append_status(out_dir, "doublet_detection", "START", "scrublet_or_fallback")
    matrix, metadata, umap = load_inputs(export_dir)
    result = run_scrublet_by_sample(matrix, metadata)
    result.to_csv(out_dir / "tables" / "scrublet_doublet_scores.tsv", sep="\t", index=False)
    write_counts(result, out_dir)
    plot_outputs(result, umap, out_dir)

    methods = ",".join(sorted(result["method"].unique()))
    n_doublets = int((result["predicted_doublet"] == "TRUE").sum())
    append_status(out_dir, "doublet_detection", "PASS", f"{n_doublets} predicted doublets; methods={methods}")
    print(f"Wrote doublet calls for {len(result)} cells; predicted_doublets={n_doublets}; methods={methods}")


if __name__ == "__main__":
    main()
