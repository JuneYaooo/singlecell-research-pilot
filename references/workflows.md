# Single-cell RNA-seq Workflow Reference

Use this reference when planning, running, or reviewing a single-cell RNA-seq project.

## Input Inspection

Check these first:

- FASTQ: sample sheet, chemistry, read structure, expected cells, genome reference, sequencing QC.
- Cell Ranger output: `filtered_feature_bc_matrix/`, `raw_feature_bc_matrix/`, `metrics_summary.csv`, `web_summary.html`, BAM if velocity is needed.
- 10x matrix: `barcodes.tsv.gz`, `features.tsv.gz` or `genes.tsv.gz`, `matrix.mtx.gz`.
- Seurat RDS: assays, reductions, metadata columns, active identities, normalization method, integration method.
- h5ad: `.obs`, `.var`, `.obsm`, `.layers`, raw counts availability, batch and condition columns.
- loom: spliced/unspliced matrices, barcode compatibility with Seurat/h5ad object.

Record whether raw counts are available. Many downstream methods need counts, not only scaled data.

## Minimal Metadata

Require or build a template with:

- `sample_id`
- `condition`
- `batch`
- `donor_id` or `patient_id`
- `tissue`
- `timepoint` if applicable
- `platform` or chemistry
- `sex`, `age`, treatment, genotype, or clinical variables when relevant

Avoid condition-level claims when each condition has only one biological replicate. In that case, frame results as exploratory.

## Standard Analysis Path

### 1. Raw Processing

Use when starting from FASTQ. Prefer the platform-standard or existing lab pipeline:

- 10x Chromium gene expression: Cell Ranger or STARsolo.
- Cost-sensitive or custom references: STARsolo, kallisto/bustools.
- Multiome or CITE-seq: use modality-aware pipelines.

Expected outputs: count matrix, web summary, alignment/QC metrics, reference version, command log.

### 2. Object Construction

Create Seurat or AnnData object from raw counts. Keep raw counts in a stable assay/layer. Add metadata before filtering so filtered-out cells can be audited if needed.

Expected outputs: raw object, metadata table, import summary.

### 3. QC and Filtering

Common metrics:

- `nFeature_RNA`
- `nCount_RNA`
- `percent.mt`
- `percent.ribo`
- hemoglobin/stress/cell-cycle scores when tissue-relevant

Set thresholds from data distributions and tissue expectations, not fixed universal cutoffs. Save pre/post filtering plots.

Recommended checks:

- Empty droplets: DropletUtils or Cell Ranger filtering.
- Doublets: DoubletFinder, scDblFinder, Scrublet.
- Ambient RNA: SoupX, DecontX, CellBender where appropriate.

Expected outputs: violin plots, scatter plots, cell count table by sample before/after filtering, filtered object.

### 4. Normalization and Feature Selection

Typical choices:

- Seurat: `NormalizeData` + `FindVariableFeatures`, or `SCTransform`.
- Scanpy: library-size normalization + log1p + highly variable genes.

Keep raw counts accessible for DE, proportion tests, and some communication methods.

### 5. Batch Assessment and Integration

Inspect UMAP/PCA colored by sample, batch, condition, and cell cycle before integration. Integrate only when technical effects dominate or when cross-sample annotation needs it.

Common methods:

- Seurat CCA/RPCA for Seurat workflows.
- Harmony for practical batch correction with readable embeddings.
- BBKNN, Scanorama, scVI for Python workflows.

Avoid integrating away true condition effects. Preserve an unintegrated object for expression-level analyses.

### 6. Clustering and Visualization

Run PCA, neighbor graph, Leiden/Louvain clustering, UMAP or t-SNE. Sweep resolution when cluster granularity affects conclusions.

Expected outputs: elbow/PC selection plots, UMAP by cluster/sample/condition, cluster size table, resolution notes.

### 7. Cell Annotation

Use a layered approach:

1. Find marker genes per cluster.
2. Compare against canonical markers and tissue biology.
3. Use automated reference tools as secondary evidence: SingleR, Azimuth, CellTypist, Symphony, scArches, or atlas mapping.
4. Validate ambiguous clusters and remove likely doublet/stress/low-quality clusters if justified.

Expected outputs: marker table, dot/violin/feature plots for marker genes, annotation table mapping cluster to cell type, uncertainty notes.

## Downstream Modules

### Differential Expression

Prefer cell type-specific comparisons. When multiple samples exist, use pseudobulk or mixed models rather than treating cells as independent biological replicates.

Tools: muscat, DESeq2/edgeR/limma on pseudobulk, MAST, Seurat `FindMarkers`, Scanpy `rank_genes_groups`.

Expected outputs: DE table with logFC, adjusted p-value, percent expression, volcano/heatmap/dot plots, methods note explaining replicate handling.

### Enrichment and Gene Set Analysis

Use DE gene lists or ranked statistics. State database version and background universe.

Tools: clusterProfiler, fgsea, gprofiler2, ReactomePA, MSigDB, Enrichr APIs when acceptable.

Expected outputs: enrichment tables, dot plots/bar plots, GSEA curves for key pathways.

### Cell Proportion Analysis

Use per-sample counts by cell type/condition. Treat simple cell-level proportions as descriptive unless replicated sample-level tests are run.

Tools: scProportionTest, scCODA, muscat, propeller, Dirichlet-multinomial or mixed models.

Expected outputs: stacked bar plots, per-sample proportion plots, test table, warning for low replicate count.

### Cell-cell Communication

Require credible annotations and enough cells per interacting type. Run separately per condition for comparison when sample design supports it.

Tools: CellChat, CellPhoneDB, NicheNet, NATMI, SingleCellSignalR.

Interpretation rules:

- Ligand-receptor predictions are hypotheses, not proof of physical interaction.
- Check whether key ligands/receptors are expressed in enough cells.
- Compare within the same preprocessing and annotation scheme across conditions.

Expected outputs: global interaction count/strength plots, pathway-level comparisons, ligand-receptor tables, bubble/chord/network plots, key pathway interpretation.

### Pseudotime and Trajectory

Use only when a continuous biological process is plausible: development, differentiation, activation, disease progression, or treatment response.

Tools: Monocle3, Slingshot, PAGA, Palantir, diffusion maps.

Root selection:

- Prefer known progenitor/naive/early time-point cells.
- If automatic root selection is used, document the criterion.
- Sensitivity-check root choice for major claims.

Expected outputs: trajectory UMAP, pseudotime distribution by condition/sample, genes varying over pseudotime, branch-specific markers, root choice note.

### RNA Velocity

Require spliced/unspliced counts from loom or compatible layers. Do not run velocity from a plain expression matrix.

Tools: velocyto, scVelo, dynamo, DeepVelo when appropriate.

Expected outputs: velocity stream/grid plots, latent time, velocity genes, PAGA if used, model mode and convergence notes.

### scDown-style Downstream Run

Use when the user wants an integrated downstream package similar to scDown:

- Input: annotated Seurat RDS or Scanpy h5ad, plus loom/spliced-unspliced files for velocity.
- Modules: annotation transfer with Symphony, cell proportion analysis with scProportionTest, communication with CellChat, pseudotime with Monocle3, RNA velocity with scVelo.
- Outputs: RDS/h5ad objects, CSV tables, high-resolution figures, per-module logs.

Confirm that conditions, cell type labels, and sample IDs are present before running.

## Report Structure

Use this structure for `report.md`:

1. Project summary
2. Input data and metadata
3. Methods and software versions
4. QC summary
5. Clustering and annotation
6. Main downstream results
7. Biological interpretation
8. Limitations and confounders
9. Reproducibility appendix with commands, parameters, and output paths

## Common Failure Checks

- Missing raw counts after importing preprocessed objects.
- Metadata columns renamed or inconsistent across samples.
- Condition confounded with batch or donor.
- Too few cells for a cell type-specific test.
- Too few biological replicates for DE or proportion claims.
- Over-integration hiding condition biology.
- Automated annotation accepted without marker validation.
- Communication analysis run on poor annotations.
- Pseudotime inferred across unrelated cell types.
- RNA velocity attempted without spliced/unspliced counts.

## Final Review Checklist

Before calling the analysis complete, verify:

- Main scripts/notebooks run from a clean session or documented environment.
- Output paths exist and contain expected figures/tables.
- Figures match the reported interpretation.
- Tables include adjusted p-values or clear descriptive labels.
- Report distinguishes validated facts from computational hypotheses.
- Package versions and commands are recorded.
