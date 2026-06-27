#!/usr/bin/env python3

from __future__ import annotations

import argparse
import os
import subprocess
from dataclasses import dataclass, field
from pathlib import Path


def default_course_root() -> Path | None:
    raw = os.environ.get("SINGLECELL_COURSE_ROOT") or os.environ.get("COURSE_ROOT")
    return Path(raw).expanduser() if raw else None


@dataclass(frozen=True)
class Step:
    name: str
    command: list[str]
    log_path: Path
    env: dict[str, str] = field(default_factory=dict)


def script(path: str) -> str:
    return str(Path(path))


def build_plan(
    *,
    course_root: Path,
    course_out: Path,
    workflow_out: Path,
    rscript: Path,
    python: Path,
    max_cells_per_sample: int,
    cibersort_max_samples: int,
    skip_course_core: bool,
    skip_cibersort: bool,
    skip_core_modules: bool,
    skip_downstream_modules: bool,
) -> list[Step]:
    logs = workflow_out / "logs" / "orchestrator"
    course_processed_rds = course_out / "objects" / "processed_multi10x_qc_harmony_cluster.rds"
    course_cibersort_dir = course_root / "24-1反卷积(CIBERSORT)"

    steps: list[Step] = []
    if not skip_course_core:
        steps.append(
            Step(
                name="course_core",
                command=[
                    str(rscript),
                    script("scripts/run_course_core.R"),
                    "--course-root",
                    str(course_root),
                    "--out",
                    str(course_out),
                    "--max-cells-per-sample",
                    str(max_cells_per_sample),
                ],
                log_path=logs / "course_core.log",
            )
        )

    if not skip_cibersort:
        steps.append(
            Step(
                name="cibersort",
                command=[
                    str(rscript),
                    script("scripts/run_course_cibersort.R"),
                    "--course-dir",
                    str(course_cibersort_dir),
                    "--seurat-rds",
                    str(course_processed_rds),
                    "--out",
                    str(course_out / "cibersort"),
                    "--max-samples",
                    str(cibersort_max_samples),
                ],
                log_path=logs / "cibersort.log",
            )
        )

    if not skip_core_modules:
        steps.extend(
            [
                Step(
                    name="object_export",
                    command=[
                        str(rscript),
                        script("scripts/run_module_export_seurat.R"),
                        "--input",
                        str(course_processed_rds),
                        "--out",
                        str(workflow_out),
                    ],
                    log_path=logs / "object_export.log",
                ),
                Step(
                    name="doublet_detection",
                    command=[
                        str(python),
                        script("scripts/run_module_scrublet.py"),
                        "--export-dir",
                        str(workflow_out / "objects" / "export"),
                        "--out",
                        str(workflow_out),
                    ],
                    log_path=logs / "doublet_detection.log",
                ),
                Step(
                    name="post_doublet_reprocess",
                    command=[
                        str(rscript),
                        script("scripts/run_module_reprocess_singlets.R"),
                        "--input",
                        str(course_processed_rds),
                        "--doublets",
                        str(workflow_out / "tables" / "scrublet_doublet_scores.tsv"),
                        "--out",
                        str(workflow_out),
                    ],
                    log_path=logs / "post_doublet_reprocess.log",
                ),
                Step(
                    name="enrichment",
                    command=[
                        str(python),
                        script("scripts/run_module_enrichment.py"),
                        "--markers",
                        str(workflow_out / "tables" / "cluster_markers_seurat.tsv"),
                        "--out",
                        str(workflow_out),
                    ],
                    log_path=logs / "enrichment.log",
                ),
                Step(
                    name="core_report",
                    command=[
                        str(python),
                        script("scripts/run_module_core_report.py"),
                        "--out",
                        str(workflow_out),
                    ],
                    log_path=logs / "core_report.log",
                ),
            ]
        )

    if not skip_downstream_modules:
        steps.extend(
            [
                Step(
                    name="export_consensus_embeddings",
                    command=[
                        str(rscript),
                        script("scripts/run_module_export_consensus_embeddings.R"),
                        "--out",
                        str(workflow_out),
                    ],
                    log_path=logs / "export_consensus_embeddings.log",
                ),
                Step(
                    name="communication",
                    command=[str(python), script("scripts/run_module_communication.py"), "--out", str(workflow_out)],
                    log_path=logs / "communication.log",
                ),
                Step(
                    name="trajectory",
                    command=[str(python), script("scripts/run_module_trajectory.py"), "--out", str(workflow_out)],
                    log_path=logs / "trajectory.log",
                ),
                Step(
                    name="cytotrace",
                    command=[str(python), script("scripts/run_module_cytotrace.py"), "--out", str(workflow_out)],
                    log_path=logs / "cytotrace.log",
                ),
                Step(
                    name="subset_extraction",
                    command=[
                        str(rscript),
                        script("scripts/run_module_subset.R"),
                        "--out",
                        str(workflow_out),
                    ],
                    log_path=logs / "subset_extraction.log",
                ),
                Step(
                    name="findmarkers_pairwise",
                    command=[
                        str(rscript),
                        script("scripts/run_module_findmarkers.R"),
                        "--out",
                        str(workflow_out),
                    ],
                    log_path=logs / "findmarkers_pairwise.log",
                ),
                Step(
                    name="gene_expression_plotting",
                    command=[
                        str(rscript),
                        script("scripts/run_module_gene_expression_plotting.R"),
                        "--out",
                        str(workflow_out),
                    ],
                    log_path=logs / "gene_expression_plotting.log",
                ),
                Step(
                    name="cnv_proxy",
                    command=[str(python), script("scripts/run_module_cnv_proxy.py"), "--out", str(workflow_out)],
                    log_path=logs / "cnv_proxy.log",
                ),
                Step(
                    name="coexpression",
                    command=[str(python), script("scripts/run_module_coexpression.py"), "--out", str(workflow_out)],
                    log_path=logs / "coexpression.log",
                ),
                Step(
                    name="celltypist_bounded",
                    command=[str(python), script("scripts/run_module_celltypist_bounded.py"), "--out", str(workflow_out)],
                    log_path=logs / "celltypist_bounded.log",
                ),
                Step(
                    name="source_traceability",
                    command=[
                        str(python),
                        script("scripts/build_source_traceability.py"),
                        "--course-root",
                        str(course_root),
                        "--workflow-dir",
                        str(workflow_out),
                    ],
                    log_path=logs / "source_traceability.log",
                ),
                Step(
                    name="file_manifest",
                    command=[str(python), script("scripts/run_module_file_manifest.py"), "--out", str(workflow_out)],
                    log_path=logs / "file_manifest.log",
                ),
                Step(
                    name="final_report",
                    command=[str(python), script("scripts/run_module_final_report.py"), "--out", str(workflow_out)],
                    log_path=logs / "final_report.log",
                ),
                Step(
                    name="final_artifact_manifest",
                    command=[
                        str(python),
                        script("scripts/run_module_file_manifest.py"),
                        "--out",
                        str(workflow_out),
                        "--final",
                    ],
                    log_path=logs / "final_artifact_manifest.log",
                ),
            ]
        )

    return steps


def run_step(step: Step) -> None:
    step.log_path.parent.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env.update(step.env)
    result = subprocess.run(step.command, text=True, capture_output=True, env=env)
    step.log_path.write_text(
        "$ " + " ".join(step.command) + "\n\n"
        + "## stdout\n"
        + result.stdout
        + "\n## stderr\n"
        + result.stderr
    )
    if result.returncode != 0:
        raise RuntimeError(f"Step {step.name} failed with exit {result.returncode}; see {step.log_path}")


def run_plan(plan: list[Step], *, dry_run: bool = False) -> None:
    for step in plan:
        print(f"=== {step.name} ===")
        print(" ".join(step.command))
        if dry_run:
            continue
        run_step(step)
        print(f"log: {step.log_path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Python orchestrator for the reference single-cell workflow.")
    parser.add_argument(
        "--course-root",
        default=default_course_root(),
        type=Path,
        help="Extracted course root. Can also be set with SINGLECELL_COURSE_ROOT.",
    )
    parser.add_argument("--course-out", default=Path("analysis/course_run"), type=Path)
    parser.add_argument("--workflow-out", default=Path("analysis/workflow_run"), type=Path)
    parser.add_argument("--rscript", default=Path(".conda/seurat-core/bin/Rscript"), type=Path)
    parser.add_argument("--python", default=Path(".conda/scverse-course/bin/python"), type=Path)
    parser.add_argument("--max-cells-per-sample", default=2500, type=int)
    parser.add_argument("--cibersort-max-samples", default=30, type=int)
    parser.add_argument("--skip-course-core", action="store_true")
    parser.add_argument("--skip-cibersort", action="store_true")
    parser.add_argument("--skip-core-modules", action="store_true")
    parser.add_argument("--skip-downstream-modules", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.course_root is None:
        raise SystemExit(
            "--course-root or SINGLECELL_COURSE_ROOT is required. "
            "Point it at the extracted course directory."
        )
    plan = build_plan(
        course_root=args.course_root,
        course_out=args.course_out,
        workflow_out=args.workflow_out,
        rscript=args.rscript,
        python=args.python,
        max_cells_per_sample=args.max_cells_per_sample,
        cibersort_max_samples=args.cibersort_max_samples,
        skip_course_core=args.skip_course_core,
        skip_cibersort=args.skip_cibersort,
        skip_core_modules=args.skip_core_modules,
        skip_downstream_modules=args.skip_downstream_modules,
    )
    run_plan(plan, dry_run=args.dry_run)
    if not args.dry_run:
        print(f"Completed reference workflow. Final report: {args.workflow_out / 'final_report.md'}")


if __name__ == "__main__":
    main()
