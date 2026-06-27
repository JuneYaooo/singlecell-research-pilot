#!/usr/bin/env python3

from __future__ import annotations

import argparse
import subprocess
from dataclasses import dataclass
from pathlib import Path

import pandas as pd


class AcceptanceError(RuntimeError):
    pass


@dataclass(frozen=True)
class CheckResult:
    check: str
    status: str
    detail: str


@dataclass(frozen=True)
class TableCheck:
    path: str
    required_columns: tuple[str, ...]
    min_data_rows: int = 1


EXPECTED_ROW_COUNTS = {
    "tables/communication_lr_scores.tsv": 641,
    "tables/trajectory_pseudotime.tsv": 2951,
    "tables/cnv_proxy_scores.tsv": 4042,
    "tables/coexpression_gene_modules.tsv": 1501,
    "tables/coexpression_module_scores.tsv": 4042,
    "tables/coexpression_module_celltype_summary.tsv": 9,
    "tables/celltypist_bounded_annotation.tsv": 2,
    "tables/final_module_status.tsv": 24,
}

REQUIRED_LOGIC_FILES = [
    "final_report.md",
    "source_traceability_report.md",
    "tables/source_traceability_matrix.tsv",
    "tables/annotated_consensus_metadata.tsv",
    "tables/annotated_consensus_umap.tsv",
    "tables/annotated_consensus_reductions.tsv",
    "figures/communication_lr_heatmap.png",
    "figures/trajectory_pseudotime_umap.png",
    "tables/cytotrace_scores.tsv",
    "tables/cytotrace_celltype_summary.tsv",
    "tables/cytotrace_method_status.tsv",
    "figures/cytotrace_umap.png",
    "objects/rds/subset_T_cells.rds",
    "tables/subset_T_cells_barcodes.tsv",
    "tables/subset_T_cells_summary.tsv",
    "tables/subset_T_cells_params.tsv",
    "tables/findmarkers_pairwise.tsv",
    "tables/findmarkers_pairwise_params.tsv",
    "figures/findmarkers_pairwise_volcano.png",
    "tables/gene_expression_requested_genes.tsv",
    "tables/gene_expression_summary.tsv",
    "figures/gene_expression_featureplot.png",
    "figures/gene_expression_vlnplot.png",
    "figures/gene_expression_dotplot.png",
    "tables/celltypist_model_status.tsv",
    "tables/enrichment_method_status.tsv",
    "figures/cnv_proxy_umap.png",
    "figures/coexpression_module_celltype_heatmap.png",
    "tables/file_manifest.tsv",
    "tables/file_manifest_summary.tsv",
    "tables/final_artifact_manifest.tsv",
    "tables/final_artifact_manifest_summary.tsv",
    "objects/rds/annotated_consensus.rds",
]

REQUIRED_COURSE_FILES = [
    "tables/run_status.tsv",
    "objects/processed_multi10x_qc_harmony_cluster.rds",
    "tables/markers_FindAllMarkers.tsv",
    "cibersort/cibersort_status.tsv",
    "cibersort/CIBERSORT-Results.tsv",
]

REQUIRED_COMPLETION_MODULES = {
    "export_consensus_embeddings",
    "cytotrace",
    "subset_extraction",
    "findmarkers_pairwise",
    "gene_expression_plotting",
    "file_manifest",
    "final_artifact_manifest",
}

FLEXIBLE_TABLE_CHECKS = [
    TableCheck(
        path="tables/annotated_consensus_umap.tsv",
        required_columns=("barcode", "umap_1", "umap_2"),
        min_data_rows=100,
    ),
    TableCheck(
        path="tables/cytotrace_method_status.tsv",
        required_columns=("module", "method", "status", "detail", "output_table"),
        min_data_rows=1,
    ),
    TableCheck(
        path="tables/cytotrace_scores.tsv",
        required_columns=(
            "barcode",
            "celltype",
            "n_detected_genes",
            "total_counts",
            "cytotrace_like_score",
            "method",
            "equivalence_level",
        ),
        min_data_rows=100,
    ),
    TableCheck(
        path="tables/subset_T_cells_barcodes.tsv",
        required_columns=("barcode", "subset_name", "subset_column", "subset_value"),
        min_data_rows=100,
    ),
    TableCheck(
        path="tables/findmarkers_pairwise.tsv",
        required_columns=("gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "comparison", "group_by"),
        min_data_rows=10,
    ),
    TableCheck(
        path="tables/gene_expression_requested_genes.tsv",
        required_columns=("gene", "present", "assay", "slot"),
        min_data_rows=1,
    ),
    TableCheck(
        path="tables/celltypist_model_status.tsv",
        required_columns=("module", "status", "model_requested", "model_path", "model_dir", "local_pkl_models", "reason", "output_table"),
        min_data_rows=1,
    ),
    TableCheck(
        path="tables/enrichment_method_status.tsv",
        required_columns=("module", "method", "status", "detail", "output_table"),
        min_data_rows=1,
    ),
    TableCheck(
        path="tables/file_manifest.tsv",
        required_columns=("path", "module", "kind", "size_bytes", "sha256"),
        min_data_rows=20,
    ),
    TableCheck(
        path="tables/final_artifact_manifest.tsv",
        required_columns=("path", "module", "kind", "size_bytes", "sha256"),
        min_data_rows=20,
    ),
]


def read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise AcceptanceError(f"missing TSV: {path}")
    return pd.read_csv(path, sep="\t")


def count_rows(path: Path) -> int:
    if not path.exists():
        raise AcceptanceError(f"missing table for row count: {path}")
    with path.open("rt") as handle:
        return sum(1 for _ in handle)


def check_status_table(path: Path, *, status_column: str, forbidden: tuple[str, ...]) -> CheckResult:
    table = read_tsv(path)
    if status_column not in table.columns:
        raise AcceptanceError(f"{path} missing status column {status_column}")
    bad = table[table[status_column].astype(str).isin(forbidden)]
    if not bad.empty:
        raise AcceptanceError(f"{path} contains forbidden statuses: {bad.to_dict(orient='records')}")
    return CheckResult(str(path), "PASS", f"{len(table)} status rows checked")


def check_file_nonempty(path: Path) -> CheckResult:
    if not path.is_file():
        raise AcceptanceError(f"required file missing: {path}")
    size = path.stat().st_size
    if size <= 0:
        raise AcceptanceError(f"required file empty: {path}")
    return CheckResult(str(path), "PASS", f"{size} bytes")


def check_row_counts(logic_out: Path) -> list[CheckResult]:
    results = []
    for rel_path, expected in EXPECTED_ROW_COUNTS.items():
        path = logic_out / rel_path
        actual = count_rows(path)
        if actual != expected:
            raise AcceptanceError(f"{path} row count {actual} != expected {expected}")
        results.append(CheckResult(str(path), "PASS", f"{actual} rows including header"))
    return results


def check_flexible_table(logic_out: Path, table_check: TableCheck) -> CheckResult:
    path = logic_out / table_check.path
    table = read_tsv(path)
    missing = set(table_check.required_columns).difference(table.columns)
    if missing:
        raise AcceptanceError(f"{path} missing required columns: {', '.join(sorted(missing))}")
    if len(table) < table_check.min_data_rows:
        raise AcceptanceError(
            f"{path} has {len(table)} data rows; expected at least {table_check.min_data_rows}"
        )
    return CheckResult(
        str(path),
        "PASS",
        f"{len(table)} data rows; columns checked: {', '.join(table_check.required_columns)}",
    )


def check_flexible_tables(logic_out: Path) -> list[CheckResult]:
    return [check_flexible_table(logic_out, table_check) for table_check in FLEXIBLE_TABLE_CHECKS]


def check_required_module_statuses(path: Path) -> CheckResult:
    table = read_tsv(path)
    required = {"module", "status", "detail"}
    missing_columns = required.difference(table.columns)
    if missing_columns:
        raise AcceptanceError(f"{path} missing columns: {', '.join(sorted(missing_columns))}")

    latest = table.drop_duplicates(subset=["module"], keep="last").set_index("module")
    missing_modules = REQUIRED_COMPLETION_MODULES.difference(latest.index.astype(str))
    if missing_modules:
        raise AcceptanceError(f"{path} missing required completion modules: {', '.join(sorted(missing_modules))}")

    allowed = {"PASS", "SKIPPED_DEPENDENCY", "SKIPPED_ASSUMPTION", "SKIPPED_REQUIRES_APPROVAL"}
    bad_rows = latest.loc[sorted(REQUIRED_COMPLETION_MODULES)]
    bad_status = bad_rows[~bad_rows["status"].astype(str).isin(allowed)]
    if not bad_status.empty:
        raise AcceptanceError(f"{path} has unacceptable completion module statuses: {bad_status.to_dict(orient='index')}")

    skipped = bad_rows[bad_rows["status"].astype(str).str.startswith("SKIPPED", na=False)]
    if not skipped.empty and skipped["detail"].fillna("").astype(str).str.strip().eq("").any():
        raise AcceptanceError(f"{path} has skipped completion modules without detail")

    return CheckResult(
        str(path),
        "PASS",
        f"completion modules present: {', '.join(sorted(REQUIRED_COMPLETION_MODULES))}",
    )


def check_file_manifest_contains_required_files(logic_out: Path) -> CheckResult:
    path = logic_out / "tables/file_manifest.tsv"
    manifest = read_tsv(path)
    if "path" not in manifest.columns:
        raise AcceptanceError(f"{path} missing path column")
    manifest_paths = set(manifest["path"].astype(str))
    expected = {
        "tables/cytotrace_scores.tsv",
        "figures/cytotrace_umap.png",
        "objects/rds/subset_T_cells.rds",
        "tables/subset_T_cells_barcodes.tsv",
        "tables/findmarkers_pairwise.tsv",
        "figures/findmarkers_pairwise_volcano.png",
        "tables/gene_expression_requested_genes.tsv",
        "figures/gene_expression_featureplot.png",
        "tables/source_traceability_matrix.tsv",
        "source_traceability_report.md",
        "tables/annotated_consensus_umap.tsv",
        "tables/cytotrace_method_status.tsv",
        "tables/celltypist_model_status.tsv",
        "tables/enrichment_method_status.tsv",
    }
    missing = expected.difference(manifest_paths)
    if missing:
        raise AcceptanceError(f"{path} missing required manifest paths: {', '.join(sorted(missing))}")
    return CheckResult(str(path), "PASS", f"{len(expected)} required paths present in manifest")


def check_traceability_matrix(path: Path) -> CheckResult:
    table = read_tsv(path)
    required = {
        "course_module",
        "source_script",
        "original_key_calls",
        "logic_preserved",
        "deviation_from_original",
        "equivalence_level",
        "acceptance_status",
        "validation_evidence",
        "python_only_feasibility",
    }
    missing = required - set(table.columns)
    if missing:
        raise AcceptanceError(f"{path} missing traceability columns: {', '.join(sorted(missing))}")
    if len(table) < 27:
        raise AcceptanceError(f"{path} has {len(table)} rows; expected at least 27")
    if table["validation_evidence"].astype(str).str.contains("MISSING|NOT_FILE", regex=True).any():
        raise AcceptanceError(f"{path} contains missing or non-file evidence")
    proxy = table[table["equivalence_level"].eq("PROXY_ONLY")]
    if not {"20_cellchat", "21_22_cnv"}.issubset(set(proxy["course_module"])):
        raise AcceptanceError("CellChat and CNV rows must be PROXY_ONLY")
    if proxy["acceptance_status"].eq("ACCEPTED").any():
        raise AcceptanceError("PROXY_ONLY rows cannot be ACCEPTED")
    return CheckResult(str(path), "PASS", f"{len(table)} traceability rows checked")


def check_celltypist_skip(path: Path) -> CheckResult:
    table = read_tsv(path)
    row = table.iloc[0].to_dict()
    expected = {
        "status": "SKIPPED_DEPENDENCY",
        "counts_cells": 4255,
        "barcodes_rows": 4255,
        "umap_rows": 4041,
        "consensus_rows": 4041,
    }
    for key, value in expected.items():
        if key not in row:
            raise AcceptanceError(f"{path} missing {key}")
        actual = row[key]
        if key != "status":
            actual = int(actual)
        if actual != value:
            raise AcceptanceError(f"{path} {key}={actual!r} != {value!r}")
    return CheckResult(str(path), "PASS", "bounded CellTypist skip metadata is consistent")


def check_rds_with_r(rscript: Path, rds_path: Path) -> CheckResult:
    if not rds_path.is_file():
        raise AcceptanceError(f"missing RDS: {rds_path}")
    code = (
        "suppressPackageStartupMessages(library(Seurat)); "
        f"obj <- readRDS('{rds_path}'); "
        "dims <- dim(obj); "
        "reds <- names(obj@reductions); "
        "if (dims[1] != 23315 || dims[2] != 4041) stop(paste('bad dims', paste(dims, collapse='x'))); "
        "if (!all(c('pca','harmony','umap') %in% reds)) stop(paste('bad reductions', paste(reds, collapse=','))); "
        "cat('RDS OK', dims[1], dims[2], paste(reds, collapse=','), '\\n')"
    )
    result = subprocess.run([str(rscript), "-e", code], text=True, capture_output=True)
    if result.returncode != 0:
        raise AcceptanceError(f"RDS check failed: {result.stderr.strip()}")
    return CheckResult(str(rds_path), "PASS", result.stdout.strip())


def write_report(results: list[CheckResult], acceptance_dir: Path) -> None:
    acceptance_dir.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame([result.__dict__ for result in results])
    frame.to_csv(acceptance_dir / "acceptance_manifest.tsv", sep="\t", index=False)
    lines = [
        "# Final Acceptance Report",
        "",
        "All checks below must pass for the current mixed R/Python course run to be considered verified.",
        "",
        "| check | status | detail |",
        "| --- | --- | --- |",
    ]
    for result in results:
        detail = result.detail.replace("|", "/")
        lines.append(f"| {result.check} | {result.status} | {detail} |")
    (acceptance_dir / "acceptance_report.md").write_text("\n".join(lines) + "\n")


def run_acceptance(course_out: Path, logic_out: Path, rscript: Path, acceptance_dir: Path) -> list[CheckResult]:
    results: list[CheckResult] = []
    results.append(check_status_table(course_out / "tables/run_status.tsv", status_column="status", forbidden=("FAIL", "ERROR")))
    results.append(check_status_table(course_out / "cibersort/cibersort_status.tsv", status_column="status", forbidden=("FAIL", "ERROR")))
    results.append(check_status_table(logic_out / "tables/final_module_status.tsv", status_column="status", forbidden=("FAIL",)))

    for rel_path in REQUIRED_COURSE_FILES:
        results.append(check_file_nonempty(course_out / rel_path))
    for rel_path in REQUIRED_LOGIC_FILES:
        results.append(check_file_nonempty(logic_out / rel_path))

    results.extend(check_row_counts(logic_out))
    results.append(check_required_module_statuses(logic_out / "tables/final_module_status.tsv"))
    results.extend(check_flexible_tables(logic_out))
    results.append(check_file_manifest_contains_required_files(logic_out))
    results.append(check_traceability_matrix(logic_out / "tables/source_traceability_matrix.tsv"))
    results.append(check_celltypist_skip(logic_out / "tables/celltypist_bounded_annotation.tsv"))
    results.append(check_rds_with_r(rscript, logic_out / "objects/rds/annotated_consensus.rds"))
    write_report(results, acceptance_dir)
    return results


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Verify final acceptance criteria for the mixed R/Python Seurat V5 run.")
    parser.add_argument("--course-out", default=Path("analysis/seurat_v5_course_run"), type=Path)
    parser.add_argument("--logic-out", default=Path("analysis/seurat_v5_logic_run"), type=Path)
    parser.add_argument("--rscript", default=Path(".conda/seurat-core/bin/Rscript"), type=Path)
    parser.add_argument("--acceptance-dir", default=None, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    acceptance_dir = args.acceptance_dir or (args.logic_out / "acceptance")
    try:
        results = run_acceptance(args.course_out, args.logic_out, args.rscript, acceptance_dir)
    except AcceptanceError as exc:
        acceptance_dir.mkdir(parents=True, exist_ok=True)
        (acceptance_dir / "acceptance_error.txt").write_text(str(exc) + "\n")
        raise SystemExit(f"ACCEPTANCE_FAIL: {exc}")
    print(f"ACCEPTANCE_PASS: {len(results)} checks passed")
    print(f"Wrote {acceptance_dir / 'acceptance_report.md'}")


if __name__ == "__main__":
    main()
