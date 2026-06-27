---
name: singlecell-research
description: Use when Codex needs to plan, run, audit, or report single-cell RNA-seq and single-cell transcriptomics analyses from FASTQ, 10x matrices, Seurat RDS, h5ad, loom, or public datasets, including QC, integration, clustering, annotation, markers, differential expression, enrichment, communication, pseudotime, CytoTRACE-style scoring, RNA velocity, CNV, co-expression, deconvolution, reproducible artifacts, and biological interpretation.
---

# Single-cell Research

## Overview

Drive a single-cell transcriptomics project from data intake to reproducible outputs and cautious biological interpretation. Treat this skill as a workflow controller: inspect inputs, choose the smallest defensible path, run or draft executable code, verify artifacts, and report what is evidence, what is statistical inference, and what is only a hypothesis.

Use the references only when needed:

| Reference | Use when |
|---|---|
| `references/workflows.md` | Planning or reviewing a general scRNA-seq workflow. |
| `references/full-workflow-contract.md` | The user asks for a complete platform-style workflow, complete module coverage, or strict acceptance criteria. |
| `references/course-adaptation.md` | The user provides the Chinese Seurat course archive or asks to adapt course-style modules. |
| `references/course-code-index.md` | You need exact source-script mapping from the course archive to an analysis step. |

## Operating Rules

- Start with provenance: dataset source, organism, tissue, disease/control design, sample count, platform/chemistry, input type, and whether raw counts are available.
- Preserve reproducibility: write commands, package versions, parameters, random seeds, output paths, and session information.
- Keep raw counts available for differential expression, communication, CNV, deconvolution, and audit.
- Prefer existing project scripts, notebooks, lockfiles, containers, and metadata over inventing a new environment.
- Protect sensitive human data. Do not upload expression matrices, metadata, or clinical files to external services unless the user explicitly approves the destination.
- Never claim biology from one visualization alone. Label results as descriptive, statistical, or hypothesis-generating.
- Ask for missing essentials only when blocked: input path, organism, metadata, biological comparison, privacy constraints, or desired output type.

## Intake

Collect or infer:

- `input_path`: files, directories, accession, object, or existing results.
- `input_type`: FASTQ, Cell Ranger output, 10x matrix, RDS, h5ad, loom, CSV/TSV matrix, or public accession.
- `organism`: human, mouse, or other.
- `sample_metadata`: sample ID, condition, batch, donor/patient, tissue, time point, platform, and replicate.
- `primary_question`: atlas construction, disease/control comparison, cell type discovery, mechanism mining, reproduction of an example, or report generation.
- `privacy_level`: public, internal, clinical/sensitive.

If metadata is missing, create a metadata template and name the required columns. Do not fabricate metadata or biological groups.

## Workflow Logic

Choose the starting point from the input:

1. FASTQ: run or plan raw processing first with Cell Ranger, STARsolo, kallisto/bustools, or an existing lab pipeline. Include read-level QC and reference version.
2. 10x or count matrix: start with object construction, metadata binding, raw-count preservation, and QC.
3. Seurat RDS or h5ad: inspect assays/layers, metadata, embeddings, clusters, normalization, integration, and annotations before rerunning upstream steps.
4. Loom or spliced/unspliced layers: consider RNA velocity only after standard QC and annotation are credible.
5. Existing marker/result tables: audit provenance and limitations before drawing conclusions.

For a full workflow with no special constraints, use this spine:

1. Data import
2. Quality control
3. Doublet and ambient RNA assessment when justified
4. Normalization and feature selection
5. Batch assessment and optional integration
6. Dimensional reduction, clustering, UMAP/t-SNE
7. Cell annotation with marker evidence
8. Marker detection and differential expression
9. Downstream modules selected by question and data suitability
10. File manifest, method-status tables, and final report

For platform-style module coverage, align outputs to `references/full-workflow-contract.md`: data import, QC, clustering, annotation, communication, pseudotime, CytoTRACE-style scoring, subset extraction, FindMarkers, FindAllMarkers, gene expression plotting, and file management.

## Module Rules

- **QC and filtering**: Use dataset-specific thresholds for genes, UMIs, mitochondrial/ribosomal percentage, empty droplets, doublets, ambient RNA, stress, and cell cycle. Save pre/post plots and cell-count tables by sample.
- **Integration and clustering**: Inspect batch structure before integration. Preserve an unintegrated object for expression-level tests. Record chosen PCs, dimensions, resolution, and integration method.
- **Annotation**: Use marker-based annotation as primary evidence. Use SingleR, Azimuth, CellTypist, Symphony, scArches, or transfer mapping as support, not blind authority. Keep uncertain labels explicit.
- **Differential expression**: Prefer cell type-specific comparisons. For multi-sample condition claims, prefer pseudobulk or models that respect biological replicates.
- **Enrichment**: Report database, organism, ID conversion, gene universe, and whether the result is over-representation, ranked GSEA, or a local fallback.
- **Cell proportions**: Use per-sample composition statistics when replicates exist. Mark cell-level-only proportions as descriptive.
- **Communication**: Require credible cell labels and enough cells per sender/receiver group. Treat ligand-receptor results as inferred signaling hypotheses.
- **Pseudotime**: Run only for plausible continua. Document subset, root choice, and sensitivity limits.
- **CytoTRACE-style scoring**: State whether official CytoTRACE/CytoTRACE2 or a gene-detection proxy ran. Do not treat proxy scores as official output.
- **RNA velocity**: Require spliced/unspliced counts. Do not run velocity from a plain expression matrix.
- **CNV**: Use only for tumor or malignant/non-malignant questions. Label expression-deviation scores as proxies unless copykat, inferCNV, or infercnvpy actually ran.
- **Co-expression and deconvolution**: Require enough cells/samples and report method-specific assumptions.

## Course Archive Adaptation

When the user supplies the course archive, treat it as a source map, not a finished pipeline. Read `references/course-adaptation.md` first, then `references/course-code-index.md` if exact script mapping is needed.

Adaptation rules:

- Replace interactive `setwd(choose.dir())` patterns with explicit input, output, and config arguments.
- Replace hard-coded metadata columns with user-provided columns.
- Choose organism-aware mitochondrial, ribosomal, annotation, GO/KEGG, communication, and CNV resources.
- Use installable current packages when original packages cannot be installed, but write method-status tables that state whether the result is method-equivalent, traceable replacement, proxy, skipped, or blocked.
- Keep source traceability for each module: original script, intended logic, adapted script/output, deviation, evidence files, and acceptance status.

## Execution Pattern

When asked to run analysis:

1. Inspect files with non-destructive commands.
2. Summarize available inputs, metadata, and blockers.
3. Choose the workflow and state why.
4. Create a run directory such as `analysis/YYYYMMDD_singlecell_<project>/`.
5. Write or reuse scripts with explicit parameters.
6. Run a small smoke test before long jobs.
7. Execute full jobs with logs.
8. Verify expected figures, tables, objects, and status files exist.
9. Draft a report with methods, results, interpretation, limitations, and reproducibility appendix.

## Deliverables

A complete run should produce:

- Reproducible scripts or notebooks.
- A parameter/config file where practical.
- `figures/` with publication-ready plots.
- `tables/` with QC counts, markers, differential expression, enrichment, cell proportions, communication, trajectory, and module-specific statistics.
- `objects/` with processed RDS/h5ad/loom outputs when storage allows.
- `logs/` with command output and package versions.
- Method-status tables for optional tools, online services, local models, fallbacks, and replacements.
- A file manifest and final artifact manifest.
- `report.md` or `report.html` with methods, results, limitations, and recommended validation.

## Verification

Before saying a workflow is complete, verify:

- The main scripts run from a documented environment.
- Expected output paths exist and non-empty artifacts are present.
- Module-status tables record pass, skip, fallback, or blocked states.
- Figures match the tables and reported interpretation.
- Raw counts and metadata needed for downstream claims are retained.
- Reports distinguish facts, statistical claims, and hypotheses.
- Package versions, commands, parameters, and caveats are recorded.

If any module used a replacement or proxy, state that in the report and do not claim method equivalence.
