#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def append_status(out_dir: Path, module: str, status: str, detail: str = "") -> None:
    status_file = out_dir / "tables" / "module_status.tsv"
    status_file.parent.mkdir(parents=True, exist_ok=True)
    if not status_file.exists():
        status_file.write_text("module\tstatus\tdetail\n")
    clean_detail = " ".join(str(detail).replace("\t", " ").splitlines())
    with status_file.open("a") as handle:
        handle.write(f"{module}\t{status}\t{clean_detail}\n")


def read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t")


def latest_status(status: pd.DataFrame) -> pd.DataFrame:
    if status.empty:
        return pd.DataFrame(columns=["module", "status", "detail"])
    rows = []
    for module, group in status.groupby("module", sort=False):
        non_start = group[group["status"] != "START"]
        rows.append((non_start if not non_start.empty else group).iloc[-1])
    return pd.DataFrame(rows).reset_index(drop=True)


def markdown_table(df: pd.DataFrame, max_rows: int | None = None) -> str:
    if df.empty:
        return "_No rows._"
    show = df.head(max_rows).copy() if max_rows else df.copy()
    show = show.fillna("")
    columns = list(show.columns)
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join(["---"] * len(columns)) + " |",
    ]
    for _, row in show.iterrows():
        values = [str(row[col]).replace("|", "/") for col in columns]
        lines.append("| " + " | ".join(values) + " |")
    if max_rows and len(df) > max_rows:
        lines.append(f"\n_Only first {max_rows} of {len(df)} rows shown._")
    return "\n".join(lines)


def file_status(path: Path) -> str:
    if not path.exists():
        return "MISSING"
    if path.is_file() and path.stat().st_size > 0:
        return f"OK ({path.stat().st_size} bytes)"
    return "EMPTY"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    out_dir = Path(args.out)
    tables = out_dir / "tables"

    append_status(out_dir, "core_report", "START", "building report")

    status_events = read_tsv(tables / "module_status.tsv")
    status_latest = latest_status(status_events)
    status_latest.to_csv(tables / "module_status_latest.tsv", sep="\t", index=False)

    equivalence = read_tsv(out_dir / "config" / "module_equivalence.tsv")
    doublet_counts = read_tsv(tables / "cell_counts_before_after_doublet.tsv")
    reprocess_summary = read_tsv(tables / "singlet_reprocess_summary.tsv")
    annotation = read_tsv(tables / "consensus_annotation.tsv")
    marker_labels = read_tsv(tables / "marker_rule_annotation.tsv")
    enrichment = read_tsv(tables / "enrichment_gprofiler.tsv")
    packages = (out_dir / "logs" / "python_packages.txt").read_text() if (out_dir / "logs" / "python_packages.txt").exists() else ""

    annotation_counts = pd.DataFrame()
    if not annotation.empty and "logic_consensus_celltype" in annotation.columns:
        annotation_counts = (
            annotation.groupby("logic_consensus_celltype", dropna=False)
            .size()
            .reset_index(name="cells")
            .sort_values("cells", ascending=False)
        )

    key_files = [
        "objects/export/counts.mtx.gz",
        "objects/export/metadata.tsv",
        "tables/scrublet_doublet_scores.tsv",
        "tables/cell_counts_before_after_doublet.tsv",
        "objects/rds/singlet_processed_clustered.rds",
        "objects/rds/annotated_consensus.rds",
        "tables/cluster_markers_seurat.tsv",
        "tables/marker_rule_annotation.tsv",
        "tables/consensus_annotation.tsv",
        "tables/enrichment_gprofiler.tsv",
        "figures/doublet_score_histogram.png",
        "figures/umap_doublet_class.png",
        "figures/umap_clusters_after_doublet.png",
        "figures/umap_marker_rule_annotation.png",
        "figures/marker_dotplot_consensus.png",
        "figures/enrichment_dotplot.png",
    ]
    file_rows = pd.DataFrame(
        [{"path": rel, "status": file_status(out_dir / rel)} for rel in key_files]
    )

    commands = [
        ".conda/seurat-core/bin/Rscript scripts/run_module_export_seurat.R --input analysis/course_run/objects/processed_multi10x_qc_harmony_cluster.rds --out analysis/workflow_run",
        "conda create -y -p .conda/scverse-course python=3.11 pip",
        ".conda/scverse-course/bin/pip install scanpy scrublet celltypist gprofiler-official matplotlib seaborn pandas scipy scikit-learn",
        ".conda/scverse-course/bin/python scripts/run_module_scrublet.py --export-dir analysis/workflow_run/objects/export --out analysis/workflow_run",
        ".conda/seurat-core/bin/Rscript scripts/run_module_reprocess_singlets.R --input analysis/course_run/objects/processed_multi10x_qc_harmony_cluster.rds --doublets analysis/workflow_run/tables/scrublet_doublet_scores.tsv --out analysis/workflow_run",
        ".conda/scverse-course/bin/python scripts/run_module_enrichment.py --markers analysis/workflow_run/tables/cluster_markers_seurat.tsv --out analysis/workflow_run --online-timeout-seconds 30",
        ".conda/scverse-course/bin/python scripts/run_module_core_report.py --out analysis/workflow_run",
    ]

    report = f"""# Single-cell Core Workflow Report

## Scope

This report summarizes the core workflow outputs used by downstream single-cell modules. It preserves the course analysis logic where a course archive is supplied, while using installable tools or documented fallbacks when exact package execution is not available.

## Input

- Verified Seurat object: `analysis/course_run/objects/processed_multi10x_qc_harmony_cluster.rds`
- Output root: `analysis/workflow_run`
- Protected R environment: `.conda/seurat-core`
- Replacement Python environment: `.conda/scverse-course`

## Commands

```bash
{chr(10).join(commands)}
```

## Python Package Availability

```text
{packages.strip()}
```

## Latest Module Status

{markdown_table(status_latest)}

## Course Module Equivalence

{markdown_table(equivalence)}

## Doublet Filtering

{markdown_table(doublet_counts)}

Scrublet was attempted through the replacement doublet module. On this processed course object, Scrublet could compute scores but could not identify an automatic threshold, so the script used a documented per-sample QC top-5-percent fallback. This preserves the course logic of removing likely high-count/high-feature doublets while keeping the run bounded and verifiable.

## Singlet Reprocessing Summary

{markdown_table(reprocess_summary)}

## Consensus Annotation Counts

{markdown_table(annotation_counts)}

## Cluster Marker-Rule Labels

{markdown_table(marker_labels)}

## Enrichment

Rows in enrichment table: {len(enrichment)}

The enrichment module first attempts online gProfiler with a bounded request timeout, then falls back to `fallback_marker_overlap` if the service is unavailable or times out. The fallback is a local, reproducible replacement artifact. It is not equivalent to a full GO/KEGG statistical enrichment test, but it preserves the course logic of interpreting cluster marker genes against biological marker categories without relying on the unavailable Bioconductor annotation stack.

## Key Output Files

{markdown_table(file_rows)}

## Downstream Modules Not Covered Here

- CellTypist automatic annotation: package is installed, but bounded model retrieval was not available; the package tried to download all 61 remote models. Marker-rule annotation is used as the primary local annotation.
- LIANA cell-cell communication: run as a downstream module after annotation review.
- infercnvpy/CNV scoring: run downstream only if epithelial/tumor-like and reference populations are credible.
- PAGA/DPT trajectory: run downstream only if the annotated cells support a plausible continuum.
- hdWGCNA/co-expression: run downstream when cell numbers and metacell design are adequate.
- Kimi/DeepSeek LLM annotation: skipped unless API credentials and upload approval are provided.

## Interpretation Limits

The core outputs are technical reproductions of the analysis logic. Marker-rule annotations are transparent but still require biological review. Fallback doublet calls and fallback enrichment are practical replacements for blocked dependencies, not numerical replicas of DoubletFinder/scDblFinder or clusterProfiler.
"""

    (out_dir / "report.md").write_text(report)
    append_status(out_dir, "core_report", "PASS", "report.md and module_status_latest.tsv written")
    print("Wrote", out_dir / "report.md")


if __name__ == "__main__":
    main()
