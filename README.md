# Single-cell Research Skill

A public Codex skill for planning, running, auditing, and reporting single-cell RNA-seq workflows. It covers the common full analysis chain: data import, quality control, doublet/ambient RNA checks, normalization, integration, clustering, annotation, marker detection, differential expression, enrichment, communication, pseudotime, CytoTRACE-style scoring, subset extraction, gene expression plotting, file manifests, and final reports.

The repository contains only reusable skill instructions, references, and runnable workflow scripts. It does not include course data, private datasets, generated analysis outputs, conda environments, model caches, binary result objects, or test scripts.

## Install

Place this directory under your Codex skills directory, or install it through your normal Codex skill/plugin workflow. The skill entrypoint is `SKILL.md`, and the public skill name is:

```text
singlecell-research
```

Minimum runtime expectations depend on which modules you run:

- Python 3.10 or newer for Python modules.
- R 4.3 or newer for Seurat-oriented modules.
- R packages commonly used by the scripts: Seurat, SeuratObject, Matrix, data.table, ggplot2, patchwork, harmony, hdf5r, R.utils, limma, e1071, preprocessCore, reshape2, and doParallel.
- Python packages commonly used by the scripts: pandas, numpy, scipy, matplotlib, seaborn, scanpy/scverse-compatible packages, scrublet, celltypist, gprofiler-official, and optional CytoTRACE-related packages.

## Repository Layout

```text
SKILL.md                         Main skill instructions
agents/openai.yaml               Codex UI metadata
references/workflows.md          General scRNA-seq workflow guidance
references/full-workflow-contract.md
                                 Complete module contract and output expectations
references/course-adaptation.md  Guidance for adapting the supplied course archive
references/course-code-index.md  Course source-script mapping
scripts/                         Reusable workflow and module runners
```

## Workflow Logic

The skill follows a decision-first workflow:

1. Inspect input type and metadata.
2. Choose the correct starting layer: FASTQ, 10x/count matrix, Seurat RDS, h5ad, loom, or existing result tables.
3. Run the core analysis spine: import, QC, normalization, batch assessment, clustering, annotation, and marker evidence.
4. Add downstream modules only when the data support them.
5. Write method-status tables when modules use optional packages, online services, local models, fallbacks, or proxy methods.
6. Produce a final report that separates reproducible evidence from biological hypotheses.

For a platform-style complete workflow, use `references/full-workflow-contract.md` as the module contract. For a general analysis, start with `references/workflows.md`. If the user provides the course archive, use `references/course-adaptation.md` and `references/course-code-index.md` to keep source logic traceable.

## Optional Reference Workflow

The original course archive is not redistributed here. If you have it locally, extract it and point the scripts at the extracted course root:

```bash
export SINGLECELL_COURSE_ROOT=/path/to/extracted/course/root
```

Run the full reference workflow when R/Python environments are already prepared:

```bash
python scripts/run_reference_workflow.py \
  --course-root "$SINGLECELL_COURSE_ROOT" \
  --course-out analysis/course_run \
  --workflow-out analysis/workflow_run \
  --rscript /path/to/Rscript \
  --python /path/to/python
```

Run downstream modules only after the core object export and annotation tables already exist:

```bash
SINGLECELL_COURSE_ROOT=/path/to/extracted/course/root \
  bash scripts/run_downstream_modules.sh analysis/workflow_run
```

For a custom project, adapt the scripts to explicit input/output paths and project metadata. Do not copy example thresholds, sample names, organism assumptions, or metadata columns blindly.

## Expected Inputs

The skill can guide workflows starting from:

- FASTQ files.
- Cell Ranger or 10x matrix directories.
- H5, Matrix Market, CSV/TSV count matrices.
- Seurat RDS objects.
- h5ad or loom objects.
- Existing marker tables, embeddings, metadata, annotation tables, or module outputs.

For reproducible analysis, provide sample metadata with at least `sample_id`, and preferably `condition`, `batch`, `donor_id`, tissue, organism, platform/chemistry, time point, and relevant clinical or treatment variables.

## Expected Outputs

A complete run should produce:

- `tables/module_status.tsv`
- `tables/final_module_status.tsv`
- `tables/source_traceability_matrix.tsv` when course-source mapping is used
- `tables/file_manifest.tsv`
- `tables/final_artifact_manifest.tsv`
- module-specific method-status tables such as `cytotrace_method_status.tsv`, `enrichment_method_status.tsv`, and `celltypist_model_status.tsv`
- reusable figures, result tables, processed objects, logs, and a final report

Completion means the artifacts exist, are reproducible from recorded commands, and carry method caveats. It does not mean every biological interpretation is experimentally proven.

## Method Boundaries

This skill is strict about method honesty:

- Local ligand-receptor scoring is not CellChat equivalence.
- Expression-deviation CNV scoring is not copykat or inferCNV equivalence.
- CytoTRACE-like gene-detection scoring is not official CytoTRACE/CytoTRACE2 unless the method-status table records a passing official run.
- CellTypist is skipped unless a local `.pkl` model is supplied or a model download is explicitly allowed.
- Online enrichment may fall back to local marker-overlap scoring; the method-status table records which path ran.
- Pseudotime is relative ordering, not measured time.
- Cross-sectional scRNA-seq does not prove causality.

## Privacy

Keep private human genomic, clinical, or unpublished data local. Do not upload expression matrices, cell metadata, or patient-derived files to external annotation services, LLM APIs, or online enrichment tools unless the data owner explicitly approves the destination and risk.

## Public Packaging

The public repository intentionally ignores generated outputs, local environments, caches, binary analysis objects, and tests:

- `analysis/`
- `.conda/`, `.venv/`, `renv/library/`
- `__pycache__/`, `.pytest_cache/`
- `tests/`, `test_*.py`
- large generated objects such as `.rds`, `.h5ad`, `.loom`, and `.mtx`

Only reusable instructions, references, and runnable analysis scripts belong in this repository.
