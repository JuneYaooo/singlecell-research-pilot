# Seurat V5 Course-Derived Workflow Reference

Use this reference when a user asks to combine this skill with the provided Seurat V5 course code, asks for a Chinese Seurat V5 full workflow, or mentions modules that match the course archive such as DoubletFinder, scDblFinder, CellChat, copykat, inferCNV, hdWGCNA, CIBERSORT, or MuSiC.

For platform-style module coverage, expected user-facing outputs, interpretation notes, and acceptance criteria, read `references/singlecell-full-workflow-user-facing.md`. That reference complements this file but does not replace source-code traceability.

The course archive is a useful implementation map, not a ready-to-run pipeline. The scripts are interactive, use `setwd(choose.dir())`, assume local example files, emit fixed output names, and often load broad package sets. Adapt their logic into project-specific scripts with explicit input paths, output directories, parameters, seeds, and session info.

For exact source-script mapping, read `references/seurat-v5-code-index.md`.

## Source Pattern

The archive covers a Seurat V5-centered R workflow:

1. Build Seurat objects from count matrices, 10X folders, multiple 10X samples, and H5 files.
2. Add mitochondrial and ribosomal QC metrics.
3. Filter cells after pre/post QC visualization.
4. Normalize, find variable features, scale, run PCA, t-SNE, UMAP, and Harmony.
5. Remove doublets with DoubletFinder or scDblFinder.
6. Cluster at multiple resolutions and choose a biologically reasonable granularity.
7. Annotate manually or with SingleR, SCINA, TransferData, scPred, or LLM-assisted marker summaries.
8. Compute markers, enrichment, pseudotime, cell-cell communication, CNV, co-expression modules, and bulk deconvolution.

Prefer `qs::qsave()` only when the environment has `qs`; otherwise use `saveRDS()` for portability.

## Input Import

Course scripts support these input styles:

- Single expression matrix: read a tabular count matrix, create a Seurat object with `CreateSeuratObject()`, then add `percent.mt` and `percent.rb`.
- Multiple matrices: create one object per dataset, add sample labels, then merge.
- Standard 10X folders: use `Read10X(data.dir = dir)` followed by `CreateSeuratObject()`.
- Non-standard 10X files: rename or decompress `barcodes`, `features`, and `matrix` files into a standard layout before import.
- Multiple 10X samples: iterate sample folders, create per-sample objects, add sample IDs, then merge.
- 10X H5 files: use `Read10X_h5()` and preserve feature names.
- Mixed formats: import each source into a Seurat object first, harmonize metadata columns, then merge.

Minimum object-construction outputs:

- Raw Seurat object before filtering.
- Metadata with `sample_id`, condition or group, batch, and donor/patient where available.
- QC metrics: `nFeature_RNA`, `nCount_RNA`, `percent.mt`, and `percent.rb`.
- Import log listing source files and object dimensions.

Guardrails:

- Do not assume human mitochondrial genes always match `^MT-`; mouse usually uses `^mt-`.
- Keep raw counts available for differential expression, CellChat, inferCNV, copykat, and deconvolution.
- Add sample metadata before merging so downstream plots and statistics can stratify by sample.

## QC and Filtering

The course QC script visualizes metrics before and after filtering:

- `VlnPlot()` for `nFeature_RNA`, `nCount_RNA`, `percent.mt`, and `percent.rb`.
- `FeatureScatter()` for count/mitochondrial, count/ribosomal, and count/feature relationships.
- Pre/post filtering PDFs.

The example thresholds are:

- `nCount_RNA <= 25000`
- `nFeature_RNA >= 0 & nFeature_RNA <= 5000`
- `percent.mt <= 25`
- `percent.rb <= 40`

Treat these as example values from the course, not universal defaults. Choose thresholds from the current dataset's distributions, tissue biology, chemistry, and expected cell types. Report how many cells and genes remain per sample after filtering.

Doublet checks:

- DoubletFinder branch: run `paramSweep()`, `summarizeSweep()`, `find.pK()`, adjust expected doublets by `modelHomotypic()`, and subset singlets.
- scDblFinder branch: convert to `SingleCellExperiment`, run `scDblFinder(samples = "sample column")`, plot `scDblFinder.class`, and subset `singlet`.

Guardrails:

- DoubletFinder metadata column names include parameter values; never hard-code a column like `DF.classifications_0.25_19_417` without discovering it first.
- Use per-sample doublet detection when samples differ in loading density.
- Keep a table of removed cells by sample and doublet class.

## Normalization, Integration, and Clustering

The course normalization branch uses:

- `NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000)`
- `FindVariableFeatures(selection.method = "vst", nfeatures = 2000)`
- `ScaleData()`
- `RunPCA()`
- `RunTSNE()` and `RunUMAP()`
- `decontX()` for ambient RNA-style count adjustment
- `RunHarmony(group.by.vars = "Type")` for batch correction in the example scripts

Adaptation rules:

- Replace `Type` with the actual batch/sample column.
- Inspect embeddings before and after integration. Do not integrate away real condition biology.
- Save an unintegrated object for expression-level testing when integration is used for visualization and clustering.
- Record chosen PCs using elbow plots, JackStraw where practical, variance explained, and biological marker coherence.

Clustering branch:

- Run `FindNeighbors(..., reduction = "harmony")` when Harmony embeddings are present.
- Sweep `FindClusters()` resolutions and visualize with `clustree`.
- Choose one resolution for downstream annotation and record the reason.
- Run `RunUMAP()` and `RunTSNE()` on the selected reduction/dimensions.

The course example includes a very low selected resolution (`0.01`) for one dataset. Do not reuse that value blindly; choose a resolution that matches sample complexity and annotation needs.

## Annotation

Manual annotation in the course uses marker dot plots, feature plots, and subclustering:

- First assign broad compartments such as immune, tumor/epithelial, stromal, endothelial, or other tissue-relevant groups.
- Re-normalize and recluster large compartments when subtype resolution is needed.
- Use `DotPlot()`, `FeaturePlot()`, `DimPlot()`, and marker tables to support labels.

Automated annotation options in the archive:

- SingleR: compare expression profiles against celldex or other reference datasets, then visualize annotation heatmaps and UMAP/t-SNE labels.
- SCINA: use marker gene lists for semi-supervised labels.
- TransferData: map query cells to a user-provided reference Seurat object.
- scPred: train/use a classifier from a labeled reference.
- LLM-assisted marker annotation: summarize markers through Kimi or DeepSeek-compatible APIs via R packages.

Guardrails:

- Automated labels are secondary evidence. Validate with canonical markers and tissue biology.
- LLM-assisted annotation must not upload private or clinical expression data without explicit user approval. Prefer sharing marker summaries, not raw matrices.
- Keep uncertain labels explicit, for example `T cell - ambiguous` or `myeloid candidate`, rather than forcing a precise subtype.

## Marker Genes and Enrichment

The marker script compares several algorithms:

- Seurat `FindAllMarkers()`
- presto `wilcoxauc()`
- COSG `cosg()`
- starTracer `searchMarker()`

Expected outputs:

- Marker CSVs per method.
- Heatmaps for top markers by cell type.
- A reconciled marker table that notes method, cluster/cell type, log fold change or AUC, adjusted p-value where available, and expression percentage.

The enrichment script uses `clusterProfiler`, `DOSE`, and `org.Hs.eg.db` for GO and KEGG:

- Run GO CC/MF/BP and KEGG for all groups and for selected cell types.
- Convert gene symbols to Entrez IDs carefully.
- State organism database and gene universe.
- Treat enrichment as functional interpretation of a gene list, not direct pathway activity proof.

Current workspace adaptation:

- `scripts/run_logic_enrichment.py` first tries online gProfiler through `gprofiler-official` with a bounded request timeout controlled by `--online-timeout-seconds`.
- If gProfiler is unavailable or times out, it writes a local marker-overlap fallback and records the method in `tables/enrichment_method_status.tsv`.
- The fallback preserves the interpretation step over marker gene lists, but it is not equivalent to `clusterProfiler` GO/KEGG statistics.

## Downstream Modules

Use downstream modules only when their assumptions are met.

### Pseudotime

Course scripts include Monocle2 and Monocle3 branches.

- Use pseudotime only for plausible continua such as differentiation, activation, disease progression, or treatment response.
- Document root-cell or root-node selection.
- Compare pseudotime distributions by sample/condition only when metadata supports it.
- Save trajectory plots, pseudotime overlays, and genes changing along trajectory.

### CellChat

Course script pattern:

- Create a CellChat object from normalized RNA data.
- Use `CellChatDB.human` in the example; switch to mouse or another database when needed.
- Optionally subset to `"Secreted Signaling"`.
- Produce interaction count/weight networks, pathway plots, heatmaps, ligand-receptor contribution plots, centrality, and communication pattern plots.

Guardrails:

- Require credible cell-type labels and enough cells per group.
- Ligand-receptor results are inferred signaling hypotheses, not proof of physical interactions.
- For condition comparison, run the same preprocessing and annotation scheme for each group.

### copykat and inferCNV

Use CNV modules for tumor/normal or malignant/non-malignant questions, not ordinary atlas annotation.

copykat:

- Use raw counts where possible.
- Inspect predicted aneuploid/diploid labels on embeddings.
- Validate predictions against known tumor markers and sample origin.

inferCNV:

- Build expression matrix, annotation file, and gene order file.
- Set reference groups deliberately. The course includes scripts with and without normal controls.
- The course comment notes 10X data often uses `cutoff = 0.1`, while one script line uses `cutoff = 1`; choose based on inferCNV documentation and input scale.
- Save CNV heatmaps, cluster labels, and score plots.

### hdWGCNA

Course script pattern:

- Select a biologically coherent subset or cell type.
- Run `SetupForWGCNA()`, `MetacellsByGroups()`, `NormalizeMetacells()`, soft-threshold selection, `ConstructNetwork()`, module eigengenes, module scores, hub genes, and module-cell associations.

Guardrails:

- Require enough cells for stable metacells and module inference.
- Record selected soft power and network parameters.
- Interpret modules as co-expression programs requiring validation.

### Deconvolution

Course scripts include CIBERSORT and MuSiC.

- Use single-cell data to build or validate cell-type signatures.
- Use bulk expression with consistent gene identifiers and normalization.
- CIBERSORT requires a signature matrix and bulk mixture matrix.
- MuSiC requires a `SingleCellExperiment` with sample and cell-type labels.

Guardrails:

- Harmonize gene symbols before deconvolution.
- Do not compare absolute fractions across platforms without caveats.
- Report uncertainty and sample-level design.

## Execution Guardrails

Before executing adapted course logic:

- Replace `setwd(choose.dir())` with explicit `input_dir`, `output_dir`, and config parameters.
- Create `figures/`, `tables/`, `objects/`, and `logs/` subdirectories.
- Avoid fixed output-name collisions by prefixing module numbers or timestamps.
- Replace hard-coded metadata columns such as `Type` with user-supplied columns.
- Discover generated metadata column names before subsetting, especially DoubletFinder columns.
- Use organism-aware mitochondrial, ribosomal, annotation, GO/KEGG, CellChat, and inferCNV resources.
- Run small smoke tests before long jobs.
- Save package versions with `sessionInfo()` or `renv`.
- Capture commands and logs for reproducibility.
- Keep private human genomic or clinical data local unless the user explicitly approves an external service.

Minimum final deliverables for a Seurat V5 run:

- Reproducible R scripts or notebooks.
- Config/parameter file.
- Processed object snapshots at major stages.
- QC, clustering, annotation, and downstream figures.
- Tables for metadata counts, markers, enrichment, communication, CNV, modules, or deconvolution as applicable.
- A report that separates descriptive results, statistical claims, and hypotheses.
