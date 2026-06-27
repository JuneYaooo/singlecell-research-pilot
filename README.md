# Single-cell Research Skill

A Codex skill for planning, running, and auditing single-cell RNA-seq workflows. It combines general scRNA-seq analysis guidance with a Seurat V5 course-code adaptation layer and a verified full-workflow contract covering import, QC, clustering, annotation, communication, pseudotime, CytoTRACE-style scoring, subset extraction, markers, gene expression plotting, file manifests, and acceptance checks.

This repository is a skill and workflow scaffold. It does not include course data, private datasets, generated analysis outputs, conda environments, or model weights.

## Install

Place this directory under your Codex skills directory, or install it through your normal Codex skill/plugin workflow. After installation, the skill entrypoint is `SKILL.md`.

Minimum runtime expectations:

- Python 3.10 or newer for the Python modules.
- R 4.3 or newer for Seurat-oriented modules.
- Seurat/SeuratObject, Matrix, data.table, ggplot2, harmony, and related R packages for course-derived R scripts.
- pandas, scipy, numpy, matplotlib, scanpy/scverse-compatible packages, and optional packages such as `celltypist`, `gprofiler-official`, and `cytotrace2-py` for Python modules.

The original Seurat V5 course archive is not redistributed here. If you want to run the course-derived example, extract the archive locally and set:

```bash
export SEURAT_V5_COURSE_ROOT=/path/to/extracted/SeuratV5/course/root
```

## Quickstart

Validate the skill structure:

```bash
python scripts/validate_seurat_v5_skill.py --root .
python /path/to/skill-creator/scripts/quick_validate.py .
```

Run the mixed R/Python course workflow when you have the extracted course data and environments prepared:

```bash
python scripts/run_full_course_logic.py \
  --course-root "$SEURAT_V5_COURSE_ROOT" \
  --course-out analysis/seurat_v5_course_run \
  --logic-out analysis/seurat_v5_logic_run \
  --rscript /path/to/Rscript \
  --python /path/to/python
```

Run Phase 2 only after Phase 1 objects and annotation tables already exist:

```bash
SEURAT_V5_COURSE_ROOT=/path/to/extracted/course \
  bash scripts/run_logic_phase2_all.sh analysis/seurat_v5_logic_run
```

For a custom project, use the references under `references/` as the workflow contract and adapt the scripts to explicit input/output paths. Do not copy course thresholds or metadata column names blindly.

## Expected Inputs

The skill can guide workflows starting from:

- FASTQ files.
- Cell Ranger or 10x matrix directories.
- H5, Matrix Market, CSV/TSV count matrices.
- Seurat RDS objects.
- h5ad or loom objects.
- Existing marker tables, metadata, embeddings, or annotation outputs.

For reproducible analysis, provide sample metadata with at least `sample_id`, and preferably `condition`, `batch`, `donor_id`, tissue, organism, and platform/chemistry.

## Verification

A complete run should produce:

- `tables/module_status.tsv` and `tables/final_module_status.tsv`
- `tables/source_traceability_matrix.tsv`
- `tables/file_manifest.tsv`
- `tables/final_artifact_manifest.tsv`
- `acceptance/acceptance_report.md`
- module-specific method-status tables such as `cytotrace_method_status.tsv`, `enrichment_method_status.tsv`, and `celltypist_model_status.tsv`

Use:

```bash
python scripts/verify_final_acceptance.py \
  --course-out analysis/seurat_v5_course_run \
  --logic-out analysis/seurat_v5_logic_run \
  --rscript /path/to/Rscript
```

Passing acceptance means delivered artifacts have reproducible paths, status rows, checks, and method caveats. It does not mean every biological interpretation is experimentally proven.

## Method Boundaries

This skill is strict about method honesty:

- Local ligand-receptor scoring is not CellChat equivalence.
- Expression-deviation CNV scoring is not copykat or inferCNV equivalence.
- CytoTRACE-like gene-detection scoring is not official CytoTRACE/CytoTRACE2 unless the method-status table records a passing official run.
- CellTypist is skipped unless a local `.pkl` model is supplied.
- Online gProfiler may fall back to local marker-overlap enrichment; the method-status table records which path ran.

Use `references/singlecell-full-workflow-user-facing.md` for the full module contract and `references/seurat-v5-workflows.md` for course-derived implementation guidance.

## Privacy

Keep private human genomic, clinical, or unpublished data local. Do not upload expression matrices, cell metadata, or patient-derived files to external annotation services, LLM APIs, or online enrichment tools unless the data owner explicitly approves the destination and risk.

## Public Packaging Notes

Generated analysis outputs, local conda environments, model caches, and binary result objects are intentionally ignored by `.gitignore`. Public releases should include the skill instructions, references, scripts, tests, and lightweight validation utilities only.
