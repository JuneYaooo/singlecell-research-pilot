# Seurat Course Code Index

Source archive:

`/path/to/course-archive.zip`

Set `SINGLECELL_COURSE_ROOT` to the extracted course root before running course-derived scripts.

ZIP member names are stored with legacy Windows naming and should be decoded by re-encoding as CP437 and decoding as GBK when they appear garbled.

Total R scripts: 35

Do not copy the archive into this skill. The archive contains Windows installers, package tarballs, public example matrices, 10X tar files, TCGA expression data, and a `scRNA.qs` object. Use those files only when the user explicitly wants to run or inspect the course examples from their original location.

## How to Use This Index

When a user asks for a Seurat workflow step, find the matching script below, then adapt its logic into a parameterized project script. Do not copy the interactive `setwd(choose.dir())` style unchanged. Confirm metadata columns, organism, input paths, and output directories before running adapted code.

## Scripts

| Stage | Exact path | Role | Key calls/packages | Notes |
|---|---|---|---|---|
| Package setup | `2.R包安装/1.基础R包安装.r` | Install CRAN/Bioconductor dependencies for the course environment. | `install.packages()`, `BiocManager::install()` | Use as dependency inventory, not as a blind install script. It can modify the user's R library heavily. |
| Package setup | `2.R包安装/2.Github安装.r` | Install GitHub and local-source dependencies. | `devtools::install_github()`, CIBERSORT, CellChat, copykat, DoubletFinder, hdWGCNA | Requires credentials/network and version choices. Prefer `renv` or documented package versions. |
| Import | `3.单矩阵读取/单矩阵读取.r` | Build a Seurat object from one expression matrix. | `CreateSeuratObject()`, `PercentageFeatureSet()`, `qs::qsave()` | Adds `percent.mt` with `^MT-` and `percent.rb` with `^RP`; adjust by organism. |
| Import | `4.多矩阵读取/多矩阵读取.r` | Build and merge Seurat objects from multiple matrix files. | `CreateSeuratObject()`, merge, QC metrics | Ensure sample IDs and metadata are added before merge. |
| Import | `5.10X数据读取(标准10X)/10X数据读取.r` | Import a standard 10X matrix folder. | `Read10X()`, `CreateSeuratObject()`, QC metrics | Assumes 10X folder layout. Record source folder and dimensions. |
| Import | `6.10X数据读取(非标准10X)/非标准10X数据读取 .R` | Import non-standard 10X data after file normalization/decompression. | `R.utils`, `CreateSeuratObject()`, QC metrics | Convert file names/layout to standard expectations before import. |
| Import | `7.多个10X数据读取/多个10X数据读取.r` | Import multiple 10X samples and merge them. | `Read10X()`, `CreateSeuratObject()`, merge | Use sample-specific metadata and avoid losing sample identity. |
| Import | `8.H5数据读取/H5读取.r` | Import 10X H5 data. | `Read10X_h5()`, `Read10X()`, `hdf5r`, `CreateSeuratObject()` | Preserve feature names and matrix type. |
| Import/merge | `9.多个不同数据合并(思路)/多个不同数据合并(思路).r` | Sketch merging heterogeneous input formats. | Seurat merge utilities, metadata harmonization | Treat as conceptual guidance; make format-specific import code explicit. |
| QC | `10.数据质控/数据质控.r` | Visualize and filter cells by QC metrics. | `VlnPlot()`, `FeatureScatter()`, `subset()` | Example thresholds are `nCount_RNA <= 25000`, `nFeature_RNA <= 5000`, `percent.mt <= 25`, `percent.rb <= 40`; re-evaluate per dataset. |
| Normalization | `11.标准化(过滤游离RNA)/标准化步骤.r` | Normalize, select variable features, scale, run PCA/t-SNE/UMAP, decontaminate, and run Harmony. | `NormalizeData()`, `FindVariableFeatures()`, `ScaleData()`, `RunPCA()`, `RunHarmony()`, `decontX()` | Replace example batch column `Type` with the user's metadata column. Preserve unintegrated counts. |
| Doublet filtering | `12-1.去除双细胞(DoubletFinder)/DoubletFinder.r` | Detect and remove doublets with DoubletFinder. | `paramSweep()`, `find.pK()`, `modelHomotypic()`, `doubletFinder()` | Discover DoubletFinder classification column dynamically; do not hard-code parameter-derived column names. |
| Doublet filtering | `12-2.去除双细胞(scDblFinder)/scDblFinder.r` | Detect and remove doublets with scDblFinder. | `SingleCellExperiment`, `scDblFinder(samples = ...)`, `DimPlot()` | Use sample-specific doublet calls when sample labels exist. |
| Post-doublet processing | `12-3.去除双细胞后常规标准化/去除双细胞后常规标准化.r` | Re-run standard normalization and embedding after doublet removal. | `NormalizeData()`, `FindVariableFeatures()`, `ScaleData()`, `RunHarmony()`, `RunUMAP()`, `RunTSNE()` | Use when doublet filtering materially changes object composition. |
| Clustering | `13.分群聚类/分群聚类.r` | Sweep clustering resolution and produce UMAP/t-SNE clusters. | `FindNeighbors()`, `FindClusters()`, `clustree`, `RunUMAP()`, `RunTSNE()` | Example selects a very low resolution (`0.01`); choose resolution by biology and marker coherence. |
| Annotation | `14.细胞手动注释/细胞手动注释.r` | Manually annotate cell types and subcluster immune cells. | `DotPlot()`, `FeaturePlot()`, `DimPlot()`, `NormalizeData()`, `FindClusters()` | Marker lists must be tissue- and organism-specific. Keep uncertain labels explicit. |
| Annotation | `15-1.细胞自动注释(singleR)/single.r` | Annotate cells with SingleR and reference data. | `SingleR`, `celldex`, heatmap, `DimPlot()` | Validate labels with markers; reference choice determines usefulness. |
| Annotation | `15-2.细胞自动注释(SCINA)/SCINA.r` | Semi-supervised annotation with marker gene lists. | `SCINA`, marker lists, `DimPlot()`, `FeaturePlot()` | Requires curated markers. Ambiguous results need manual review. |
| Annotation | `15-3.细胞自动注释(调用Kimi)/调用Kimi进行自动注释.r` | LLM-assisted annotation through Kimi-compatible API. | `openai`, `xxdAIcelltype`, `FindAllMarkers()` | Do not upload private data without explicit approval. Prefer marker summaries over raw matrices. |
| Annotation | `15-4.细胞自动注释(调用DeepSeek)/调用DeepSeek进行自动注释.r` | LLM-assisted annotation through DeepSeek-compatible API. | `openai`, `xxdAIcelltype`, `FindAllMarkers()` | Same privacy guardrail as Kimi. Treat output as annotation suggestion. |
| Reference mapping | `16-1.自有模型细胞自动注释(TransferData)/TransferData.R` | Transfer labels from a user-provided Seurat reference. | `FindTransferAnchors()`, `TransferData()` | Requires a compatible labeled reference object. Check feature overlap. |
| Reference mapping | `16-2.自有模型细胞自动注释(scPred)/scPred.r` | Annotate with a trained scPred model. | `scPred`, `getFeatureSpace()`, `trainModel()`, `scPredict()` | Requires training labels and careful train/query compatibility. |
| Markers | `17.多种算法计算marker基因/多种算法计算marker基因.R` | Compare marker genes across several algorithms. | `FindAllMarkers()`, `presto::wilcoxauc()`, `COSG::cosg()`, `starTracer::searchMarker()` | Use multiple methods to prioritize robust markers; keep method-specific statistics. |
| Enrichment | `18.GO_KEGG/GO_KEGG.r` | Run GO and KEGG enrichment from marker genes. | `clusterProfiler`, `DOSE`, `org.Hs.eg.db`, `enrichGO()`, `enrichKEGG()` | Change organism database for non-human data. Report gene universe and ID conversion losses. |
| Pseudotime | `19-1.拟时序分析Monocle2/拟时序分析Monocle2.r` | Monocle2 trajectory analysis. | `monocle`, `FindAllMarkers()`, trajectory plots | Use only for plausible continua and document root choice. |
| Pseudotime | `19-2.拟时序分析Monocle3/Monocle3.r` | Monocle3 trajectory analysis and gene trends. | `monocle3`, `as.cell_data_set()`, `learn_graph()`, `order_cells()` | Requires careful cell subset and root-node/root-cell choice. |
| Communication | `20.CellChat/CellChat.r` | Cell-cell communication inference and visualization. | `CellChat`, `CellChatDB.human`, `subsetDB()`, `computeCommunProb()`, network plots | Switch database by organism. Interpret ligand-receptor predictions as hypotheses. |
| CNV | `21.copykat/copykat.r` | Infer aneuploid/diploid status with copykat. | `copykat`, Seurat visualization | Best suited to tumor datasets with plausible CNV signal. Validate with biology. |
| CNV | `22-1.infercnv/infercnv.r` | Run inferCNV without explicit normal-control script branch. | `CreateInfercnvObject()`, `infercnv::run()`, `plot_cnv()` | Set reference groups deliberately; choose cutoff based on data scale. |
| CNV | `22-2.infercnv(引入正常对照)/infercnv(引入正常对照).r` | Run inferCNV using normal-control cells. | Seurat preprocessing, `CreateInfercnvObject()`, `infercnv::run()` | Requires credible normal reference annotations. |
| CNV | `22-3.infercnv(计算打分)/infercnv(计算打分).r` | Score and visualize inferCNV-derived malignant/normal patterns. | `ComplexHeatmap`, `circlize`, inferCNV object parsing, boxplots | Use as scoring-pattern reference; confirm object paths and group labels. |
| Co-expression | `23.hdWGCNA/hdWGCNA.r` | Build hdWGCNA modules and visualize module programs. | `hdWGCNA`, `SetupForWGCNA()`, `MetacellsByGroups()`, `ConstructNetwork()`, hub genes | Requires enough cells/metacells in the selected population. Record soft power. |
| Deconvolution helper | `24-1反卷积(CIBERSORT)/xxdimmune.R` | CIBERSORT helper implementation. | `CIBERSORT()`, `e1071`, `preprocessCore`, `limma` | Treat as bundled algorithm code. Check licensing and provenance before reuse. |
| Deconvolution | `24-1反卷积(CIBERSORT)/反卷积(CIBERSORT).R` | Run CIBERSORT-style bulk deconvolution. | `CIBERSORT`, signature matrix, bulk expression matrix, boxplots | Requires harmonized genes and a valid signature matrix. |
| Deconvolution | `24-2反卷积(MusiC)/反卷积(MusiC).r` | Run MuSiC deconvolution using single-cell reference and bulk data. | `MuSiC`, `SingleCellExperiment`, `music_prop()` | Requires sample and cell-type labels in the single-cell reference. |

## Module-to-Script Lookup

- Need import code: use scripts 3 through 9.
- Need QC thresholds and plots: use script 10.
- Need ambient RNA/decontX and Harmony: use script 11.
- Need doublet filtering: use scripts 12-1 or 12-2.
- Need clustering resolution sweep: use script 13.
- Need annotation: use scripts 14 through 16-2.
- Need marker and enrichment tables: use scripts 17 and 18.
- Need trajectory: use scripts 19-1 or 19-2.
- Need CellChat: use script 20.
- Need malignant CNV calls: use scripts 21 and 22-1 through 22-3.
- Need co-expression modules: use script 23.
- Need bulk deconvolution: use scripts 24-1 and 24-2.
