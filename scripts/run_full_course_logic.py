#!/usr/bin/env python3

from __future__ import annotations

import argparse
import os
import subprocess
from dataclasses import dataclass, field
from pathlib import Path


def default_course_root() -> Path | None:
    raw = os.environ.get("SEURAT_V5_COURSE_ROOT") or os.environ.get("COURSE_ROOT")
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
    logic_out: Path,
    rscript: Path,
    python: Path,
    max_cells_per_sample: int,
    cibersort_max_samples: int,
    skip_course_core: bool,
    skip_cibersort: bool,
    skip_phase1: bool,
    skip_phase2: bool,
) -> list[Step]:
    logs = logic_out / "logs" / "orchestrator"
    course_processed_rds = course_out / "objects" / "processed_multi10x_qc_harmony_cluster.rds"
    course_cibersort_dir = course_root / "24-1反卷积(CIBERSORT)"

    steps: list[Step] = []
    if not skip_course_core:
        steps.append(
            Step(
                name="course_core",
                command=[
                    str(rscript),
                    script("scripts/run_seurat_v5_course_core.R"),
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
                    script("scripts/run_seurat_v5_course_cibersort.R"),
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

    if not skip_phase1:
        steps.extend(
            [
                Step(
                    name="logic_export",
                    command=[
                        str(rscript),
                        script("scripts/run_logic_export_seurat.R"),
                        "--input",
                        str(course_processed_rds),
                        "--out",
                        str(logic_out),
                    ],
                    log_path=logs / "logic_export.log",
                ),
                Step(
                    name="doublet_detection",
                    command=[
                        str(python),
                        script("scripts/run_logic_scrublet.py"),
                        "--export-dir",
                        str(logic_out / "objects" / "export"),
                        "--out",
                        str(logic_out),
                    ],
                    log_path=logs / "doublet_detection.log",
                ),
                Step(
                    name="post_doublet_reprocess",
                    command=[
                        str(rscript),
                        script("scripts/run_logic_reprocess_singlets.R"),
                        "--input",
                        str(course_processed_rds),
                        "--doublets",
                        str(logic_out / "tables" / "scrublet_doublet_scores.tsv"),
                        "--out",
                        str(logic_out),
                    ],
                    log_path=logs / "post_doublet_reprocess.log",
                ),
                Step(
                    name="enrichment",
                    command=[
                        str(python),
                        script("scripts/run_logic_enrichment.py"),
                        "--markers",
                        str(logic_out / "tables" / "cluster_markers_seurat.tsv"),
                        "--out",
                        str(logic_out),
                    ],
                    log_path=logs / "enrichment.log",
                ),
                Step(
                    name="phase1_report",
                    command=[
                        str(python),
                        script("scripts/run_logic_report.py"),
                        "--out",
                        str(logic_out),
                    ],
                    log_path=logs / "phase1_report.log",
                ),
            ]
        )

    if not skip_phase2:
        steps.extend(
            [
                Step(
                    name="export_consensus_embeddings",
                    command=[
                        str(rscript),
                        script("scripts/run_logic_export_consensus_embeddings.R"),
                        "--out",
                        str(logic_out),
                    ],
                    log_path=logs / "export_consensus_embeddings.log",
                ),
                Step(
                    name="communication",
                    command=[str(python), script("scripts/run_logic_communication.py"), "--out", str(logic_out)],
                    log_path=logs / "communication.log",
                ),
                Step(
                    name="trajectory",
                    command=[str(python), script("scripts/run_logic_trajectory.py"), "--out", str(logic_out)],
                    log_path=logs / "trajectory.log",
                ),
                Step(
                    name="cytotrace",
                    command=[str(python), script("scripts/run_logic_cytotrace.py"), "--out", str(logic_out)],
                    log_path=logs / "cytotrace.log",
                ),
                Step(
                    name="subset_extraction",
                    command=[
                        str(rscript),
                        script("scripts/run_logic_subset.R"),
                        "--out",
                        str(logic_out),
                    ],
                    log_path=logs / "subset_extraction.log",
                ),
                Step(
                    name="findmarkers_pairwise",
                    command=[
                        str(rscript),
                        script("scripts/run_logic_findmarkers.R"),
                        "--out",
                        str(logic_out),
                    ],
                    log_path=logs / "findmarkers_pairwise.log",
                ),
                Step(
                    name="gene_expression_plotting",
                    command=[
                        str(rscript),
                        script("scripts/run_logic_gene_expression_plotting.R"),
                        "--out",
                        str(logic_out),
                    ],
                    log_path=logs / "gene_expression_plotting.log",
                ),
                Step(
                    name="cnv_proxy",
                    command=[str(python), script("scripts/run_logic_cnv_proxy.py"), "--out", str(logic_out)],
                    log_path=logs / "cnv_proxy.log",
                ),
                Step(
                    name="coexpression",
                    command=[str(python), script("scripts/run_logic_coexpression.py"), "--out", str(logic_out)],
                    log_path=logs / "coexpression.log",
                ),
                Step(
                    name="celltypist_bounded",
                    command=[str(python), script("scripts/run_logic_celltypist_bounded.py"), "--out", str(logic_out)],
                    log_path=logs / "celltypist_bounded.log",
                ),
                Step(
                    name="source_traceability",
                    command=[
                        str(python),
                        script("scripts/build_source_traceability.py"),
                        "--course-root",
                        str(course_root),
                        "--logic-dir",
                        str(logic_out),
                    ],
                    log_path=logs / "source_traceability.log",
                ),
                Step(
                    name="file_manifest",
                    command=[str(python), script("scripts/run_logic_file_manifest.py"), "--out", str(logic_out)],
                    log_path=logs / "file_manifest.log",
                ),
                Step(
                    name="final_report",
                    command=[str(python), script("scripts/run_logic_final_report.py"), "--out", str(logic_out)],
                    log_path=logs / "final_report.log",
                ),
                Step(
                    name="final_artifact_manifest",
                    command=[
                        str(python),
                        script("scripts/run_logic_file_manifest.py"),
                        "--out",
                        str(logic_out),
                        "--final",
                    ],
                    log_path=logs / "final_artifact_manifest.log",
                ),
            ]
        )

    steps.append(
        Step(
            name="final_acceptance",
            command=[
                str(python),
                script("scripts/verify_final_acceptance.py"),
                "--course-out",
                str(course_out),
                "--logic-out",
                str(logic_out),
                "--rscript",
                str(rscript),
            ],
            log_path=logs / "final_acceptance.log",
        )
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
    parser = argparse.ArgumentParser(description="Python orchestrator for the mixed R/Python Seurat V5 course workflow.")
    parser.add_argument(
        "--course-root",
        default=default_course_root(),
        type=Path,
        help="Extracted Seurat V5 course root. Can also be set with SEURAT_V5_COURSE_ROOT.",
    )
    parser.add_argument("--course-out", default=Path("analysis/seurat_v5_course_run"), type=Path)
    parser.add_argument("--logic-out", default=Path("analysis/seurat_v5_logic_run"), type=Path)
    parser.add_argument("--rscript", default=Path(".conda/seurat-core/bin/Rscript"), type=Path)
    parser.add_argument("--python", default=Path(".conda/scverse-course/bin/python"), type=Path)
    parser.add_argument("--max-cells-per-sample", default=2500, type=int)
    parser.add_argument("--cibersort-max-samples", default=30, type=int)
    parser.add_argument("--skip-course-core", action="store_true")
    parser.add_argument("--skip-cibersort", action="store_true")
    parser.add_argument("--skip-phase1", action="store_true")
    parser.add_argument("--skip-phase2", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.course_root is None:
        raise SystemExit(
            "--course-root or SEURAT_V5_COURSE_ROOT is required. "
            "Point it at the extracted Seurat V5 course directory."
        )
    plan = build_plan(
        course_root=args.course_root,
        course_out=args.course_out,
        logic_out=args.logic_out,
        rscript=args.rscript,
        python=args.python,
        max_cells_per_sample=args.max_cells_per_sample,
        cibersort_max_samples=args.cibersort_max_samples,
        skip_course_core=args.skip_course_core,
        skip_cibersort=args.skip_cibersort,
        skip_phase1=args.skip_phase1,
        skip_phase2=args.skip_phase2,
    )
    run_plan(plan, dry_run=args.dry_run)
    if not args.dry_run:
        print(f"Completed mixed R/Python course workflow. Final report: {args.logic_out / 'final_report.md'}")


if __name__ == "__main__":
    main()
