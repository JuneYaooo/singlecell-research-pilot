# Single-cell Full Workflow User-facing Reference

Use this reference when the user asks for a complete, reusable, no-code-style or platform-style single-cell transcriptomics workflow. This file converts the user-provided Zero Coding platform article into a workflow contract: what each module is for, what inputs it needs, what outputs it should produce, how it maps to the Seurat course material, and how to validate it.

This is a product/workflow reference, not method evidence. Use `references/course-code-index.md` and `references/course-adaptation.md` for course-source traceability. Use `analysis/workflow_run/tables/source_traceability_matrix.tsv`, final status tables, manifests, and reports for implementation evidence.

## Design Boundary

The full skill has three complementary layers:

| layer | source | role |
|---|---|---|
| User workflow contract | This reference, based on the user-provided platform article | Defines modules, expected outputs, interpretation notes, and validation criteria. |
| Method traceability | Seurat course scripts and code index | Anchors each module to original commands where available. |
| Executable evidence | Local R/Python scripts, status tables, manifests, and reports | Proves what currently runs in this environment and labels replacements honestly. |

Do not treat a user-facing module list as proof that a method ran. A module is complete only when its inputs, command path, logs, outputs, and artifact checks are present.

## Recommended Module Order

For a standard 10x or count-matrix scRNA-seq dataset with sample metadata, run:

1. Data Import
2. Quality Control
3. Clustering
4. Cell Annotation
5. Cell-cell Communication
6. Pseudotime
7. CytoTRACE
8. Subset Extraction
9. FindMarkers
10. FindAllMarkers
11. Gene Expression Plotting
12. File Management

The first four modules are the analysis spine. Downstream modules must use the object, metadata, cell labels, and quality decisions produced by those core steps.

## Verified Workspace Artifacts

The bundled reference implementation exports the final annotated Seurat object into reusable audit tables before downstream Python modules run:

- `tables/annotated_consensus_metadata.tsv`
- `tables/annotated_consensus_umap.tsv`
- `tables/annotated_consensus_reductions.tsv`

Method-status files are required for modules that may use optional packages, online services, local models, or traceable replacements:

- `tables/cytotrace_method_status.tsv` records whether `cytotrace2_py` was available, whether official CytoTRACE2 execution was attempted, and which score table was delivered.
- `tables/enrichment_method_status.tsv` records whether online gProfiler succeeded or local marker-overlap fallback was used. The enrichment command supports `--online-timeout-seconds` so external API calls cannot block the workflow indefinitely.
- `tables/celltypist_model_status.tsv` records the requested CellTypist model path, local model discovery, and skip reason when no `.pkl` model is available.

Final audit files must exist before claiming a full run is complete:

- `tables/final_module_status.tsv`
- `tables/final_artifact_manifest.tsv`
- `tables/final_artifact_manifest_summary.tsv`
- `final_report.md`

## Module Contract

### Data Import

Purpose: create the first analysis object from FASTQ-derived matrices, Cell Ranger output, 10x folders, H5 files, Seurat RDS, h5ad, loom, or tabular count matrices.

Prerequisites:

- Input path and file type.
- Organism and expected mitochondrial gene prefix.
- Sample metadata with at least `sample_id`; add `condition`, `batch`, and `donor_id` when available.

Course mapping:

- Scripts 3 through 9 in `references/course-code-index.md`.
- Key logic: `Read10X()`, `Read10X_h5()`, `CreateSeuratObject()`, merge, metadata harmonization, `PercentageFeatureSet()`.

Required outputs:

- Raw object before filtering.
- Import summary with cells, genes, samples, input paths, and organism.
- Metadata table with sample labels.
- Initial QC figures for `nFeature_RNA`, `nCount_RNA`, and `percent.mt`.

Validation:

- Raw counts are preserved.
- Cell barcodes are unique after merge.
- Sample labels are not missing.
- Import dimensions match source files or expected Cell Ranger summary.

Interpretation notes:

- Use this step to judge whether data are analyzable, not to draw biology conclusions.
- High mitochondrial percentage or extreme UMI counts should guide QC thresholds.

### Quality Control

Purpose: remove low-quality cells, likely empty droplets, high-mitochondrial cells, abnormal UMI/gene outliers, ambient RNA where justified, and likely doublets.

Prerequisites:

- Imported object with raw counts and metadata.
- Clear QC thresholds or data-driven threshold selection rules.

Course mapping:

- QC script 10.
- Normalization/decontX script 11 for ambient RNA-style logic.
- Doublet scripts 12-1, 12-2, and 12-3.

Required outputs:

- Pre-QC and post-QC violin plots.
- Pre-QC and post-QC scatter/correlation plots.
- Cell count table by sample before and after filtering.
- Filtered object.
- Doublet table if doublet detection was run.

Validation:

- Removed cell counts are reported by sample and reason.
- Thresholds are recorded and not copied blindly from the course.
- Post-QC object has raw counts, metadata, and enough cells for downstream analysis.

Interpretation notes:

- QC thresholds are dataset-specific.
- Overly strict filtering can remove real rare cell states.
- Doublet and ambient RNA calls are quality-control evidence, not cell-type labels.

### Clustering

Purpose: identify transcriptionally similar cell groups and generate low-dimensional visualizations.

Prerequisites:

- QC-filtered object.
- Normalization method and batch/integration decision.
- PC/dimension selection rule.

Course mapping:

- Script 11 for normalization, variable features, scaling, PCA, UMAP/t-SNE, Harmony.
- Script 13 for neighbor graph, clustering, and resolution sweep.

Required outputs:

- Highly variable gene plot or table.
- Elbow/PC selection figure.
- PCA plot colored by sample, condition, and batch where available.
- Batch-corrected embedding plot if integration is used.
- UMAP and t-SNE by cluster.
- Cluster size table.
- Selected resolution note.

Validation:

- Embeddings and cluster labels exist in the object.
- Chosen PCs and resolution are recorded.
- Batch correction is justified and not used blindly.
- UMAP/t-SNE plots are consistent with cluster tables.

Interpretation notes:

- Clustering defines the structure used by annotation and downstream modules.
- Too many clusters can split one cell type into technical subgroups; too few can hide real subtypes.

### Cell Annotation

Purpose: translate numeric clusters into biologically interpretable cell types or states.

Prerequisites:

- Clustered object.
- Marker table or marker gene list.
- Organism, tissue context, and disease/control design.

Course mapping:

- Manual annotation script 14.
- Automated annotation scripts 15-1 through 16-2.
- Marker script 17.

Required outputs:

- Cluster-to-cell-type annotation table.
- Marker table supporting each label.
- UMAP and t-SNE colored by cell type.
- Cell-type plots faceted or grouped by condition/sample.
- Cell-type composition plots.
- Annotation uncertainty notes.

Validation:

- Each major label has marker evidence.
- Automated labels are checked against markers and tissue biology.
- Ambiguous clusters remain labeled as ambiguous instead of being forced.
- Sample and condition composition tables are available.

Interpretation notes:

- Automated annotation is supporting evidence, not final authority.
- Cell composition differences require sample-level design for strong claims.

### Cell-cell Communication

Purpose: infer candidate ligand-receptor interactions between annotated cell types.

Prerequisites:

- Credible cell-type labels.
- Enough cells per sender and receiver group.
- Organism-appropriate ligand-receptor database.

Course mapping:

- CellChat script 20.
- The bundled reference implementation uses curated ligand-receptor scoring and is a proxy, not CellChat equivalence.

Required outputs:

- Global network count/weight figures.
- Per-cell-type or pathway network figures.
- Bubble/chord/heatmap plots where supported.
- Ligand-receptor score table.
- Pathway-level summary table when method supports it.

Validation:

- Sender and receiver labels match annotation.
- Ligands and receptors are present in expression data.
- Results are labeled as inferred signaling hypotheses.
- Proxy methods are not marked as CellChat-equivalent.

Interpretation notes:

- Communication results do not prove physical interaction.
- Compare conditions only when preprocessing and annotation are consistent across groups.

### Pseudotime

Purpose: model a plausible continuous cell-state process such as differentiation, activation, tumor progression, or treatment response.

Prerequisites:

- A biologically coherent subset or cell-type group.
- Root-cell/root-node strategy.
- Metadata for condition/sample comparisons.

Course mapping:

- Monocle2 script 19-1.
- Monocle3 script 19-2.
- Current verified replacement uses Scanpy diffusion pseudotime/PAGA and must be labeled as a traceable replacement, not Monocle output.

Required outputs:

- Pseudotime trajectory plot.
- State/branch plot when method supports it.
- Trajectory colored by cell type.
- Trajectory colored by condition/sample.
- Pseudotime table per cell.
- Root choice note.

Validation:

- Cell subset is documented.
- Root choice is documented.
- Pseudotime values exist for all included cells.
- Claims are limited to relative ordering unless time-course design supports more.

Interpretation notes:

- Pseudotime is relative order, not measured time.
- Do not infer trajectories across unrelated cell types.

### CytoTRACE

Purpose: score relative differentiation potential and compare it with cluster, annotation, and pseudotime patterns.

Prerequisites:

- Gene expression matrix with adequate gene detection.
- Cell metadata and embeddings.
- Appropriate species/gene identifier compatibility.

Course mapping:

- The provided Seurat course index does not currently include a dedicated CytoTRACE script.
- Treat this as a user-facing workflow requirement from the platform article and implement separately if needed.

Required outputs:

- CytoTRACE score table per cell.
- UMAP/t-SNE colored by CytoTRACE score.
- Score summaries by cluster, cell type, and condition.
- Comparison note against pseudotime if both are run.

Validation:

- Input gene IDs are compatible with the CytoTRACE implementation.
- Scores cover the expected cell set.
- Module status clearly says whether original CytoTRACE, CytoTRACE2, or another potential-scoring replacement was used.

Bundled implementation status:

- Executable as `scripts/run_module_cytotrace.py`.
- The current implementation writes `tables/cytotrace_scores.tsv`, `tables/cytotrace_celltype_summary.tsv`, `tables/cytotrace_method_status.tsv`, and `figures/cytotrace_umap.png`.
- The bundled implementation can detect `cytotrace2_py`, but official CytoTRACE2 execution is disabled by default because that API can require dense file materialization and model downloads. The delivered default score remains a CytoTRACE-like gene-detection potential score, not official CytoTRACE equivalence.
- If official CytoTRACE2 is explicitly attempted later, `tables/cytotrace_method_status.tsv` must show `method`, `status`, `official_package_available`, and `official_execution` rather than relying on narrative notes.

Interpretation notes:

- Higher score usually indicates greater inferred developmental potential.
- CytoTRACE and pseudotime answer related but different questions.

### Subset Extraction

Purpose: extract one cluster, cell type, sample group, condition, or custom metadata subset into a reusable object for focused downstream analysis.

Prerequisites:

- Annotated object.
- Explicit subset rule.
- Decision on whether to re-normalize/recluster the subset.

Course mapping:

- Subclustering logic appears in manual annotation script 14 and downstream modules that select populations.
- No standalone productized subset module is currently validated.

Required outputs:

- Subset object.
- Subset summary table with cells by sample, condition, cluster, and cell type.
- Subset rule log.

Validation:

- Subset rule is reproducible.
- Extracted barcodes match the rule.
- Parent object and subset object paths are recorded.

Bundled implementation status:

- Executable as `scripts/run_module_subset.R`.
- The current default extracts `logic_consensus_celltype == T_cells` from the annotated Seurat object and writes `objects/rds/subset_T_cells.rds`, barcode membership, parameters, and summary tables.
- The module is parameterized with `--subset-column`, `--subset-values`, and `--name` for other datasets.

Interpretation notes:

- Subset extraction is an enabling step, not a biological result by itself.

### FindMarkers

Purpose: compare two specified groups for differential expression.

Prerequisites:

- Explicit group A and group B.
- Cell type, cluster, condition, or metadata comparison definition.
- Replicate-aware design if claiming condition-level biology.

Course mapping:

- Seurat marker logic from script 17.
- Pairwise `FindMarkers()` is part of Seurat logic and is exposed as a standalone executable module in the bundled implementation.

Required outputs:

- Pairwise differential expression table.
- Parameters: assay/layer, test method, logFC threshold, minimum expression percentage, covariates if used.
- Volcano, heatmap, dot plot, or top-gene table as appropriate.

Validation:

- Both groups contain enough cells.
- The comparison is not confounded with one sample unless clearly exploratory.
- Adjusted p-values and effect sizes are present when method supports them.

Bundled implementation status:

- Executable as `scripts/run_module_findmarkers.R`.
- The current default compares `T_cells` versus `B_cells` by `logic_consensus_celltype` and writes `tables/findmarkers_pairwise.tsv`, parameters, and a volcano plot.
- The module is parameterized with `--group-by`, `--ident-1`, `--ident-2`, and optional subset arguments.

Interpretation notes:

- For multi-sample disease/control claims, prefer pseudobulk or replicate-aware models over treating all cells as independent.

### FindAllMarkers

Purpose: find marker genes for every cluster or every annotated cell type.

Prerequisites:

- Cluster or cell-type identities.
- Normalized expression and raw counts where required.

Course mapping:

- Script 17: `FindAllMarkers()`, `presto::wilcoxauc()`, `COSG::cosg()`, `starTracer::searchMarker()`.
- The bundled reference workflow includes Seurat `FindAllMarkers()` on the singlet object.

Required outputs:

- Marker table for all clusters or cell types.
- Top-marker heatmap or dot plot.
- Reconciled marker summary if multiple algorithms are used.

Validation:

- Marker rows are nonempty.
- Group labels match object metadata.
- Method, thresholds, and adjusted p-value rules are recorded.

Interpretation notes:

- FindAllMarkers identifies representative genes, not necessarily condition differential genes.

### Gene Expression Plotting

Purpose: visualize selected gene expression on embeddings or by group for marker checking, candidate-gene display, and result validation.

Prerequisites:

- Processed object with embeddings.
- Gene list and grouping variable.
- Confirmed gene symbols/IDs.

Course mapping:

- Manual annotation and marker scripts use `FeaturePlot()`, `VlnPlot()`, `DotPlot()`, `DimPlot()`, and heatmaps.

Required outputs:

- FeaturePlot on UMAP and/or t-SNE.
- Violin plot by cluster, cell type, condition, or sample.
- DotPlot for marker panels.
- Missing-gene report.

Validation:

- Requested genes exist or are reported as missing.
- Plots use the intended assay/layer.
- Figure captions state whether values are raw, normalized, scaled, or imputed.

Bundled implementation status:

- Executable as `scripts/run_module_gene_expression_plotting.R`.
- The current default plots `CD3D,MS4A1,EPCAM,LYZ,COL1A1,VWF` from the annotated Seurat object and writes requested-gene status, expression summaries, FeaturePlot, VlnPlot, and DotPlot outputs.
- The module is parameterized with `--genes`, `--group-by`, `--assay`, and `--slot`.

Interpretation notes:

- A marker plot validates annotation only when expression pattern matches the expected cell population and metadata context.

### File Management

Purpose: make results reproducible, downloadable, and auditable.

Prerequisites:

- Run directory.
- Stable subdirectories for figures, tables, objects, logs, and reports.

Course mapping:

- The course scripts write many fixed-name files; adapted runs must parameterize output directories and avoid collisions.
- The bundled implementation has reports, module status tables, source traceability, final manifests, and a SHA256 file manifest module.

Required outputs:

- `figures/`
- `tables/`
- `objects/`
- `logs/`
- `report.md` or `report.html`
- Package/session info.
- Manifest with file paths, sizes, checksums where practical, and module ownership.
- Final manifest and final status snapshots: `tables/final_artifact_manifest.tsv`, `tables/final_artifact_manifest_summary.tsv`, and `tables/final_module_status.tsv`.

Validation:

- Expected files exist and are nonempty.
- Module status has no hidden failures.
- Artifact checks can be rerun from the recorded commands and manifests.
- Skipped modules record why they were skipped.
- `final_artifact_manifest.tsv` includes final reports and manifest files, so the downloadable artifact set is auditable by checksum.

Interpretation notes:

- File management is part of scientific reproducibility, not a cosmetic download feature.

## Validation Matrix

Use this matrix before saying a workflow is complete.

| module | minimum completion signal | bundled implementation status |
|---|---|---|
| Data Import | Object created, metadata added, QC metrics present | Verified through course core and logic export. |
| Quality Control | Pre/post QC outputs and filtered object | Verified for course example; doublet replacement has caveat. |
| Clustering | PCA/integration/UMAP clusters present | Verified with Seurat/Harmony object. |
| Cell Annotation | Cell-type labels and marker evidence | Verified as broad marker-rule consensus. |
| Cell-cell Communication | LR tables and network/heatmap outputs | Proxy only; not CellChat-equivalent. |
| Pseudotime | Pseudotime table and trajectory plot | Traceable replacement; not Monocle-equivalent. |
| CytoTRACE | Score table, embedding plot, and method-status table | Executable `cytotrace2_py`-aware CytoTRACE-like replacement; not official CytoTRACE-equivalent unless status says official run passed. |
| Subset Extraction | Standalone subset object and rule log | Executable default T-cell subset with parameterized subset rules. |
| FindMarkers | Pairwise DE module | Executable Seurat `FindMarkers()` module. |
| FindAllMarkers | All-group marker table | Verified with Seurat FindAllMarkers. |
| Gene Expression Plotting | Arbitrary-gene plotting module | Executable Seurat FeaturePlot, VlnPlot, and DotPlot module. |
| File Management | Manifest, final manifest, reports, logs, final status | Executable SHA256 manifest, `final_artifact_manifest.tsv`, `final_module_status.tsv`, reports, logs, and artifact checks. |

## Replacement Policy

Use these labels consistently:

| label | meaning | claim allowed |
|---|---|---|
| COURSE_CONSISTENT | Same course method or same package family with parameterized execution | Can be accepted with normal validation. |
| METHOD_EQUIVALENT | Different wrapper but same statistical method and input/output contract | Can be accepted with caveat. |
| TRACEABLE_REPLACEMENT | Different method preserving the analysis question | Can be accepted as replacement, not exact reproduction. |
| PROXY_ONLY | Heuristic or bounded local proxy for a blocked method | Do not claim method equivalence. |
| SKIPPED_DEPENDENCY | Tool or model unavailable | Report as skipped, not failed biology. |
| SKIPPED_REQUIRES_APPROVAL | External API or private-data upload would be required | Ask before running. |
| BLOCKED_ORIGINAL | Original method cannot run in current environment | Keep blocked status and provide reason. |

Hard rule: proxy modules cannot satisfy original-package equivalence. For example, curated ligand-receptor scoring is not CellChat, and expression-deviation CNV scoring is not copykat or inferCNV.

## Completion Standard

A completed full-workflow skill should provide:

- A module status table covering all modules in the recommended order.
- A source traceability matrix for every course-derived or replaced method.
- A report that separates descriptive results, statistical claims, and hypotheses.
- Re-runnable commands or orchestration scripts.
- Method-status tables for optional or replacement modules, including `cytotrace_method_status.tsv`, `enrichment_method_status.tsv`, and `celltypist_model_status.tsv`.
- A final report and manifests that cover required objects, exported embeddings, figures, tables, row counts, checksums, and skipped-module reasons.

When the user asks for "100% verifiable", answer in terms of this completion standard: 100% verifiable means every delivered artifact has a reproducible source, a status, and a check. It does not mean every biological hypothesis is true or every replacement is equivalent to the original R package.
