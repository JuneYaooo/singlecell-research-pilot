---
name: singlecell-research-pilot
description: End-to-end single-cell RNA-seq research workflow pilot for Codex. Use when the user asks to analyze scRNA-seq or single-cell transcriptomics data, run Seurat/Seurat V5/Scanpy/scDown-style pipelines, combine Chinese Seurat V5 course code archives, process 10x matrices/FASTQ/RDS/h5ad/loom files, perform QC, doublet or ambient RNA filtering, normalization, integration, clustering, annotation, marker/DE/enrichment analysis, CellChat/CellPhoneDB communication, Monocle/Slingshot/PAGA pseudotime, scVelo/RNA velocity, copykat/inferCNV, hdWGCNA, CIBERSORT/MuSiC deconvolution, or produce reproducible figures, tables, scripts, and reports.
---

# Single-cell Research Pilot

## Overview

Drive a single-cell RNA-seq project from data intake to reproducible outputs and cautious biological interpretation. Treat the skill as a research workflow controller: inspect inputs, choose the smallest defensible analysis path, run or draft executable code, validate outputs, and produce a report that separates evidence from hypotheses.

For detailed module choices and output expectations, read `references/workflows.md` before planning or running a non-trivial analysis. For user-facing full-workflow modules, output contracts, and no-code/platform-style acceptance criteria, read `references/singlecell-full-workflow-user-facing.md`. For Seurat V5 course-derived implementation guidance, read `references/seurat-v5-workflows.md`. For exact mapping from the provided course archive to source scripts, read `references/seurat-v5-code-index.md`.

## Operating Rules

- Start with data provenance: record dataset source, organism, tissue, disease/control design, sample count, chemistry/platform, and whether the input is FASTQ, 10x matrix, Seurat RDS, h5ad, loom, or precomputed figures/tables.
- Prefer existing project scripts, notebooks, containers, and lockfiles over inventing a new environment.
- Never claim a biological conclusion from one visualization alone. State whether the result is descriptive, statistical, or a hypothesis for validation.
- Preserve reproducibility: save commands, package versions, parameters, random seeds, output paths, and session info.
- Protect sensitive human data. Do not upload private genomic or clinical data to external services unless the user explicitly approves and the destination is appropriate.
- Ask for missing essentials only when they block the workflow: biological comparison, sample metadata, organism/reference genome, input location, or desired output type.

## Intake Checklist

Collect or infer these before analysis:

- `input_path`: files or directories to inspect.
- `input_type`: FASTQ, Cell Ranger output, 10x matrix, RDS, h5ad, loom, CSV/TSV matrix, or public accession.
- `organism`: human, mouse, or other.
- `sample_metadata`: sample ID, condition, batch, donor/patient, time point, tissue, replicate.
- `primary_question`: atlas construction, disease/control comparison, cell type discovery, mechanism mining, or figure/report generation.
- `privacy_level`: public, internal, clinical/sensitive.

If metadata is absent, create a metadata template and tell the user exactly which columns are required.

## Workflow Decision Tree

1. If the input is FASTQ, plan raw processing first with Cell Ranger, STARsolo, kallisto/bustools, or an existing pipeline. Do not skip read-level QC.
2. If the input is a 10x matrix or count matrix, start at object construction and cell/gene QC.
3. If the input is Seurat RDS or h5ad, inspect existing metadata, embeddings, clusters, normalization, and annotations before rerunning upstream steps.
4. If the user wants "full process" and has no special constraints, run the standard path: QC -> doublet/ambient RNA checks -> normalization -> HVG -> PCA -> batch assessment/integration -> clustering -> UMAP/t-SNE -> annotation -> marker validation -> downstream modules -> report.
5. If the user wants a no-code/platform-style full workflow, align deliverables to `references/singlecell-full-workflow-user-facing.md`: data import, QC, clustering, annotation, communication, pseudotime, CytoTRACE, subset extraction, FindMarkers, FindAllMarkers, gene expression plotting, and file management.
6. If the user cites scDown-style downstream analysis, prioritize annotated Seurat/h5ad input and run or reproduce: cell proportion differences, cell-cell communication, pseudotime, and RNA velocity where required inputs exist.

## Analysis Modules

Use modules only when the data and study design justify them:

- **QC and filtering**: mitochondrial/ribosomal percentage, detected genes, UMI counts, empty droplets, doublets, ambient RNA, cell cycle if relevant.
- **Integration and clustering**: Seurat CCA/RPCA, Harmony, BBKNN, Scanorama, scVI, Leiden/Louvain resolution sweep, batch mixing diagnostics.
- **Annotation**: marker-based annotation as the primary evidence; use SingleR, Azimuth, CellTypist, Symphony, or reference mapping as support, not blind authority.
- **Differential expression**: compare within cell types when possible; account for sample/replicate structure; prefer pseudobulk for multi-sample disease/control claims.
- **Functional enrichment**: GO, KEGG, Reactome, Hallmark, GSEA, MSigDB or organism-specific databases; report background gene universe.
- **Cell proportion**: use per-sample composition tests when biological replicates exist; mark permutation-only or cell-level tests as exploratory if sample-level replication is weak.
- **Cell-cell communication**: use CellChat, CellPhoneDB, NicheNet, NATMI, or SingleCellSignalR; interpret ligand-receptor predictions as inferred signaling hypotheses.
- **Pseudotime/trajectory**: use Monocle3, Slingshot, PAGA, or Palantir only for biologically plausible continua; document root choice.
- **RNA velocity**: use scVelo/velocyto only when spliced/unspliced counts or loom files exist; report model mode and limitations.

## Execution Pattern

When asked to run an analysis:

1. Inspect files with non-destructive commands.
2. Summarize available inputs and missing metadata.
3. Choose the workflow and state why.
4. Create a run directory such as `analysis/YYYYMMDD_singlecell_<project>/`.
5. Write or reuse scripts/notebooks with explicit parameters.
6. Run the smallest smoke test before long jobs.
7. Execute full jobs, capturing logs.
8. Validate expected outputs exist and inspect summaries.
9. Draft the result report with figures, tables, methods, limitations, and next steps.

## Deliverables

For a completed analysis, produce:

- Reproducible scripts or notebooks.
- A parameter/config file where practical.
- `figures/` with publication-ready PNG/PDF/SVG outputs.
- `tables/` with marker genes, DE results, enrichment results, cell counts/proportions, and module-specific statistics.
- `objects/` with processed RDS/h5ad/loom outputs when storage allows.
- `logs/` with command output and package versions.
- `report.md` or `report.html` containing methods, results, interpretation, limitations, and recommended validation.
- Method-status tables and a final artifact manifest when optional tools, fallbacks, or replacement modules are used.

## Reporting Standards

Write results in a cautious research style:

- Use "suggests", "is consistent with", or "candidate mechanism" for inferred results.
- Do not overstate causality from cross-sectional scRNA-seq.
- Distinguish sample-level findings from cell-level descriptive patterns.
- Mention confounders: donor imbalance, batch effects, low cell counts, tissue dissociation bias, annotation uncertainty, ambient RNA, and doublets.
- Include enough methods detail for another analyst to reproduce the result.

## References

- Read `references/workflows.md` for module-specific guidance, recommended outputs, and common failure checks.
- Read `references/singlecell-full-workflow-user-facing.md` when the user asks for a complete no-code-style workflow, platform-like module coverage, output interpretation notes, or "100% verifiable" acceptance criteria.
- Read `references/seurat-v5-workflows.md` when the user mentions Seurat V5, the provided course archive, Chinese course code, or Seurat-centered modules such as DoubletFinder, scDblFinder, CellChat, copykat, inferCNV, hdWGCNA, CIBERSORT, or MuSiC.
- Read `references/seurat-v5-code-index.md` when mapping a requested workflow step to the exact course script path or checking which archive script covers a module.
