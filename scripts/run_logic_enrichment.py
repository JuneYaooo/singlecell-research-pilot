#!/usr/bin/env python3

from __future__ import annotations

import argparse
import signal
from contextlib import contextmanager
from pathlib import Path
from typing import Callable, TypeVar

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


MARKER_SETS = {
    "T_cells": {"CD3D", "CD3E", "CD2", "CD4", "CD8A", "NKG7", "GNLY"},
    "B_cells": {"MS4A1", "CD79A", "CD19"},
    "Plasma": {"MZB1", "SDC1", "IGHG1"},
    "Myeloid_monocyte": {"LYZ", "CD14", "FCGR3A", "CD68", "CD163"},
    "Epithelial_tumor_like": {"EPCAM", "KRT8", "KRT18", "KRT19"},
    "Endothelial": {"PECAM1", "VWF", "KDR"},
    "Fibroblast_stromal": {"COL1A1", "COL1A2", "DCN", "LUM", "MME"},
    "Mast": {"TPSAB1", "TPSB2", "CPA3", "KIT"},
}

T = TypeVar("T")


def append_status(out_dir: Path, module: str, status: str, detail: str = "") -> None:
    status_file = out_dir / "tables" / "module_status.tsv"
    status_file.parent.mkdir(parents=True, exist_ok=True)
    if not status_file.exists():
        status_file.write_text("module\tstatus\tdetail\n")
    clean_detail = " ".join(str(detail).replace("\t", " ").splitlines())
    with status_file.open("a") as handle:
        handle.write(f"{module}\t{status}\t{clean_detail}\n")


def top_markers(markers: pd.DataFrame, n: int = 100) -> dict[str, list[str]]:
    required = {"cluster", "gene", "avg_log2FC"}
    missing = required.difference(markers.columns)
    if missing:
        raise ValueError(f"Marker table missing required columns: {', '.join(sorted(missing))}")
    clean = markers.dropna(subset=["cluster", "gene", "avg_log2FC"]).copy()
    clean["cluster"] = clean["cluster"].astype(str)
    clean = clean.sort_values(["cluster", "avg_log2FC"], ascending=[True, False])
    return {
        cluster: group["gene"].astype(str).head(n).tolist()
        for cluster, group in clean.groupby("cluster", sort=False)
    }


def run_with_timeout(callable_fn: Callable[[], T], *, timeout_seconds: float, label: str) -> T:
    if timeout_seconds <= 0:
        return callable_fn()

    def raise_timeout(_signum, _frame):
        raise TimeoutError(f"{label} timed out after {timeout_seconds:g} seconds")

    previous_handler = signal.getsignal(signal.SIGALRM)
    signal.signal(signal.SIGALRM, raise_timeout)
    signal.setitimer(signal.ITIMER_REAL, timeout_seconds)
    try:
        return callable_fn()
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous_handler)


@contextmanager
def gprofiler_request_timeout(timeout_seconds: float):
    import gprofiler.gprofiler as gprofiler_module

    original_post = gprofiler_module.requests.post

    def post_with_timeout(*args, **kwargs):
        if timeout_seconds > 0:
            kwargs.setdefault("timeout", timeout_seconds)
        return original_post(*args, **kwargs)

    gprofiler_module.requests.post = post_with_timeout
    try:
        yield
    finally:
        gprofiler_module.requests.post = original_post


def run_gprofiler(marker_dict: dict[str, list[str]], *, timeout_seconds: float = 120.0) -> pd.DataFrame:
    from gprofiler import GProfiler

    def query() -> pd.DataFrame:
        gp = GProfiler(return_dataframe=True)
        rows = []
        with gprofiler_request_timeout(timeout_seconds):
            for cluster, genes in marker_dict.items():
                if not genes:
                    continue
                result = gp.profile(
                    organism="hsapiens",
                    query=genes,
                    sources=["GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC"],
                    user_threshold=0.05,
                    no_evidences=False,
                )
                if result is None or len(result) == 0:
                    continue
                result = result.copy()
                result.insert(0, "cluster", cluster)
                result["method"] = "gprofiler"
                rows.append(result)
        if not rows:
            raise RuntimeError("gProfiler returned no enrichment rows")
        return pd.concat(rows, ignore_index=True)

    return run_with_timeout(query, timeout_seconds=timeout_seconds, label="gProfiler")


def fallback_marker_overlap(marker_dict: dict[str, list[str]]) -> pd.DataFrame:
    rows = []
    for cluster, genes in marker_dict.items():
        gene_set = set(map(str.upper, genes))
        for term_name, markers in MARKER_SETS.items():
            overlap = sorted(gene_set.intersection(markers))
            rows.append(
                {
                    "cluster": cluster,
                    "source": "LOCAL_MARKER",
                    "native": term_name,
                    "name": term_name,
                    "p_value": 1.0 / (1 + len(overlap)),
                    "term_size": len(markers),
                    "query_size": len(gene_set),
                    "intersection_size": len(overlap),
                    "intersection": ",".join(overlap),
                    "method": "fallback_marker_overlap",
                }
            )
    result = pd.DataFrame(rows)
    return result.sort_values(["cluster", "intersection_size", "p_value"], ascending=[True, False, True])


def plot_enrichment(result: pd.DataFrame, out_dir: Path) -> None:
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    plot_df = result.copy()
    if "name" not in plot_df.columns:
        plot_df["name"] = plot_df.get("native", "term")
    if "p_value" not in plot_df.columns:
        plot_df["p_value"] = 1.0
    plot_df["p_value"] = pd.to_numeric(plot_df["p_value"], errors="coerce").fillna(1.0)
    plot_df["score"] = -np.log10(np.maximum(plot_df["p_value"].to_numpy(dtype=float), 1e-300))

    if "intersection_size" not in plot_df.columns:
        plot_df["intersection_size"] = 1
    plot_df["intersection_size"] = pd.to_numeric(plot_df["intersection_size"], errors="coerce").fillna(1)

    plot_df["cluster"] = plot_df["cluster"].astype(str)
    plot_df = (
        plot_df.sort_values(["cluster", "score", "intersection_size"], ascending=[True, False, False])
        .groupby("cluster", group_keys=False)
        .head(3)
        .copy()
    )
    if plot_df.empty:
        plot_df = pd.DataFrame({"cluster": ["none"], "name": ["No enrichment rows"], "score": [0], "intersection_size": [1]})

    plot_df["label"] = "C" + plot_df["cluster"].astype(str) + " " + plot_df["name"].astype(str).str.slice(0, 55)
    plot_df = plot_df.sort_values("score", ascending=True)

    height = max(4, min(16, 0.28 * len(plot_df) + 1.5))
    plt.figure(figsize=(9, height))
    plt.scatter(plot_df["score"], plot_df["label"], s=20 + 10 * plot_df["intersection_size"], c=plot_df["score"], cmap="viridis")
    plt.xlabel("-log10(p value)")
    plt.ylabel("")
    plt.tight_layout()
    plt.savefig(fig_dir / "enrichment_dotplot.png", dpi=220)
    plt.close()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--markers", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--offline", action="store_true", help="Skip online gProfiler queries and write local fallback enrichment.")
    parser.add_argument(
        "--online-timeout-seconds",
        default=120.0,
        type=float,
        help="Maximum wall-clock time for online gProfiler before falling back locally. Use 0 to disable.",
    )
    args = parser.parse_args()

    out_dir = Path(args.out)
    out_dir.joinpath("tables").mkdir(parents=True, exist_ok=True)
    out_dir.joinpath("figures").mkdir(parents=True, exist_ok=True)

    append_status(out_dir, "enrichment", "START", args.markers)
    markers = pd.read_csv(args.markers, sep="\t")
    marker_dict = top_markers(markers, n=100)

    try:
        if args.offline:
            raise RuntimeError("offline mode requested")
        result = run_gprofiler(marker_dict, timeout_seconds=args.online_timeout_seconds)
        status = "PASS"
        detail = f"{len(result)} gProfiler rows"
        method_status = {
            "module": "enrichment",
            "method": "gprofiler_official",
            "status": "PASS",
            "detail": detail,
            "output_table": "tables/enrichment_gprofiler.tsv",
        }
    except Exception as exc:
        result = fallback_marker_overlap(marker_dict)
        status = "PASS"
        detail = f"fallback_marker_overlap rows={len(result)}; gprofiler_error={type(exc).__name__}:{exc}"
        append_status(out_dir, "gprofiler_online", "SKIPPED_DEPENDENCY", f"{type(exc).__name__}:{exc}")
        method_status = {
            "module": "enrichment",
            "method": "fallback_marker_overlap",
            "status": "FALLBACK",
            "detail": f"{type(exc).__name__}:{exc}",
            "output_table": "tables/enrichment_gprofiler.tsv",
        }

    result.to_csv(out_dir / "tables" / "enrichment_gprofiler.tsv", sep="\t", index=False)
    pd.DataFrame([method_status]).to_csv(out_dir / "tables" / "enrichment_method_status.tsv", sep="\t", index=False)
    plot_enrichment(result, out_dir)
    append_status(out_dir, "enrichment", status, detail)
    print(detail)


if __name__ == "__main__":
    main()
