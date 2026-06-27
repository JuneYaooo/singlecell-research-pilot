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


def set_latest_status(latest: pd.DataFrame, module: str, status: str, detail: str) -> pd.DataFrame:
    row = {"module": module, "status": status, "detail": detail}
    if latest.empty or "module" not in latest.columns:
        return pd.DataFrame([row], columns=["module", "status", "detail"])

    updated = latest.copy()
    mask = updated["module"] == module
    if mask.any():
        updated.loc[mask, ["status", "detail"]] = [status, detail]
    else:
        updated = pd.concat([updated, pd.DataFrame([row])], ignore_index=True)
    return updated


def markdown_table(df: pd.DataFrame, max_rows: int | None = None) -> str:
    if df.empty:
        return "_No rows._"
    show = df.head(max_rows).fillna("")
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


def file_status(out_dir: Path, rel: str, expected_skip: str = "") -> dict[str, str]:
    path = out_dir / rel
    if not path.exists():
        status = expected_skip or "MISSING"
    elif path.is_file() and path.stat().st_size > 0:
        status = f"OK ({path.stat().st_size} bytes)"
    else:
        status = "EMPTY"
    return {"path": rel, "status": status}


def summarize_optional_table(path: Path, label: str) -> str:
    df = read_tsv(path)
    if df.empty:
        return f"- {label}: no table written."
    return f"- {label}: {len(df)} rows, columns: {', '.join(df.columns[:8])}"


def write_report(out_dir: Path) -> None:
    tables = out_dir / "tables"

    status_events = read_tsv(tables / "module_status.tsv")
    latest = set_latest_status(
        latest_status(status_events),
        "final_report",
        "PASS",
        "final_report.md and final_module_status.tsv written",
    )
    latest.to_csv(tables / "final_module_status.tsv", sep="\t", index=False)
    latest_status_by_module = dict(zip(latest["module"], latest["status"]))

    equivalence = read_tsv(out_dir / "config" / "module_equivalence.tsv")
    traceability = read_tsv(tables / "source_traceability_matrix.tsv")
    traceability_counts = pd.DataFrame()
    if not traceability.empty and {"equivalence_level", "acceptance_status"}.issubset(traceability.columns):
        traceability_counts = (
            traceability.groupby(["equivalence_level", "acceptance_status"], dropna=False)
            .size()
            .reset_index(name="modules")
            .sort_values(["equivalence_level", "acceptance_status"])
        )
    annotation = read_tsv(tables / "consensus_annotation.tsv")
    annotation_counts = pd.DataFrame()
    if not annotation.empty and "logic_consensus_celltype" in annotation.columns:
        annotation_counts = (
            annotation.groupby("logic_consensus_celltype", dropna=False)
            .size()
            .reset_index(name="cells")
            .sort_values("cells", ascending=False)
        )

    celltypist_skip = ""
    if latest_status_by_module.get("celltypist_bounded") == "SKIPPED_DEPENDENCY":
        celltypist_skip = "EXPECTED_SKIP (celltypist_bounded SKIPPED_DEPENDENCY)"
    key_files = [
        ("report.md", ""),
        ("tables/module_status_latest.tsv", ""),
        ("tables/final_module_status.tsv", ""),
        ("tables/source_traceability_matrix.tsv", ""),
        ("source_traceability_report.md", ""),
        ("tables/annotated_consensus_metadata.tsv", ""),
        ("tables/annotated_consensus_umap.tsv", ""),
        ("tables/annotated_consensus_reductions.tsv", ""),
        ("tables/communication_lr_scores.tsv", ""),
        ("figures/communication_lr_heatmap.png", ""),
        ("tables/trajectory_pseudotime.tsv", ""),
        ("figures/trajectory_pseudotime_umap.png", ""),
        ("tables/cytotrace_scores.tsv", ""),
        ("tables/cytotrace_celltype_summary.tsv", ""),
        ("tables/cytotrace_method_status.tsv", ""),
        ("figures/cytotrace_umap.png", ""),
        ("objects/rds/subset_T_cells.rds", ""),
        ("tables/subset_T_cells_barcodes.tsv", ""),
        ("tables/subset_T_cells_summary.tsv", ""),
        ("tables/subset_T_cells_params.tsv", ""),
        ("tables/findmarkers_pairwise.tsv", ""),
        ("tables/findmarkers_pairwise_params.tsv", ""),
        ("figures/findmarkers_pairwise_volcano.png", ""),
        ("tables/gene_expression_requested_genes.tsv", ""),
        ("tables/gene_expression_summary.tsv", ""),
        ("figures/gene_expression_featureplot.png", ""),
        ("figures/gene_expression_vlnplot.png", ""),
        ("figures/gene_expression_dotplot.png", ""),
        ("tables/celltypist_model_status.tsv", ""),
        ("tables/enrichment_method_status.tsv", ""),
        ("tables/cnv_proxy_scores.tsv", ""),
        ("figures/cnv_proxy_umap.png", ""),
        ("tables/coexpression_gene_modules.tsv", ""),
        ("tables/coexpression_module_scores.tsv", ""),
        ("tables/coexpression_module_celltype_summary.tsv", ""),
        ("figures/coexpression_module_celltype_heatmap.png", ""),
        ("tables/celltypist_bounded_annotation.tsv", ""),
        ("figures/umap_celltypist_bounded.png", celltypist_skip),
        ("objects/rds/annotated_consensus.rds", ""),
        ("tables/cluster_markers_seurat.tsv", ""),
        ("tables/enrichment_gprofiler.tsv", ""),
        ("figures/enrichment_dotplot.png", ""),
        ("tables/file_manifest.tsv", ""),
        ("tables/file_manifest_summary.tsv", ""),
        ("tables/final_artifact_manifest.tsv", ""),
        ("tables/final_artifact_manifest_summary.tsv", ""),
    ]
    file_rows = pd.DataFrame([file_status(out_dir, rel, expected_skip) for rel, expected_skip in key_files])

    optional_summaries = "\n".join(
        [
            summarize_optional_table(tables / "communication_lr_scores.tsv", "Communication LR scores"),
            summarize_optional_table(tables / "trajectory_pseudotime.tsv", "Trajectory pseudotime"),
            summarize_optional_table(tables / "cytotrace_scores.tsv", "CytoTRACE-like scores"),
            summarize_optional_table(tables / "cytotrace_method_status.tsv", "CytoTRACE method status"),
            summarize_optional_table(tables / "subset_T_cells_barcodes.tsv", "Subset extraction barcodes"),
            summarize_optional_table(tables / "findmarkers_pairwise.tsv", "Pairwise FindMarkers"),
            summarize_optional_table(tables / "gene_expression_requested_genes.tsv", "Gene expression requested genes"),
            summarize_optional_table(tables / "celltypist_model_status.tsv", "CellTypist model status"),
            summarize_optional_table(tables / "enrichment_method_status.tsv", "Enrichment method status"),
            summarize_optional_table(tables / "cnv_proxy_scores.tsv", "CNV proxy"),
            summarize_optional_table(tables / "coexpression_gene_modules.tsv", "Co-expression gene modules"),
            summarize_optional_table(tables / "coexpression_module_scores.tsv", "Co-expression module scores"),
            summarize_optional_table(tables / "celltypist_bounded_annotation.tsv", "CellTypist bounded annotation"),
            summarize_optional_table(tables / "file_manifest.tsv", "File manifest"),
            summarize_optional_table(tables / "final_artifact_manifest.tsv", "Final artifact manifest"),
        ]
    )

    report = f"""# Seurat V5 Course Logic Final Report

## Scope

This report consolidates the verified Phase 1 workflow and the Phase 2 replacement modules. The objective is to reproduce the course logic with tools that can run in this environment, not to claim numerical identity with the original package outputs.

## Output Root

`{out_dir}`

## Final Module Status

{markdown_table(latest)}

## Course Module Equivalence

{markdown_table(equivalence)}

## Source Traceability Summary

{markdown_table(traceability_counts)}

## Consensus Annotation Counts

{markdown_table(annotation_counts)}

## Phase 2 Table Summaries

{optional_summaries}

## Key Files

{markdown_table(file_rows)}

## Important Caveats

- Doublet detection used the documented per-sample QC top-5-percent fallback because Scrublet could not determine an automatic threshold on the processed course object.
- Enrichment tries online gProfiler with a bounded request timeout; this run used local marker-overlap fallback because the online request timed out or was unavailable.
- `celltypist_annotation` is a legacy Phase 1 attempt; the authoritative bounded CellTypist module is `celltypist_bounded`, which only uses local `.pkl` models and skips otherwise.
- Source traceability is now audited in `source_traceability_report.md`; modules labeled `PROXY_ONLY` are not accepted as original-package-equivalent results.
- CNV output is an exploratory expression-deviation proxy, not inferCNV or copykat.
- Communication output is a local ligand-receptor scoring replacement, not CellChat.
- Trajectory output is a diffusion pseudotime replacement, not Monocle.
- CytoTRACE output is a gene-count-rank differentiation-potential replacement. `cytotrace2_py` can be detected, but official execution is disabled by default because it may require dense input materialization and model downloads; do not claim official CytoTRACE/CytoTRACE2 equivalence unless `cytotrace_method_status.tsv` records a passing official run.
- Subset extraction, pairwise FindMarkers, and gene expression plotting are now executable user-facing modules backed by the annotated Seurat object.
- Co-expression output is a bounded local correlation/module-scoring replacement, not hdWGCNA.
"""

    (out_dir / "final_report.md").write_text(report)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    out_dir = Path(args.out)
    append_status(out_dir, "final_report", "START", "building final report")
    try:
        write_report(out_dir)
        append_status(out_dir, "final_report", "PASS", "final_report.md and final_module_status.tsv written")
        write_report(out_dir)
    except Exception as exc:
        append_status(out_dir, "final_report", "FAIL", f"{type(exc).__name__}: {exc}")
        raise
    print(f"Wrote {out_dir / 'final_report.md'}")


if __name__ == "__main__":
    main()
