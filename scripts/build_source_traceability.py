#!/usr/bin/env python3

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import pandas as pd


COURSE_CONSISTENT = "COURSE_CONSISTENT"
METHOD_EQUIVALENT = "METHOD_EQUIVALENT"
TRACEABLE_REPLACEMENT = "TRACEABLE_REPLACEMENT"
PROXY_ONLY = "PROXY_ONLY"
SKIPPED_DEPENDENCY = "SKIPPED_DEPENDENCY"
SKIPPED_REQUIRES_APPROVAL = "SKIPPED_REQUIRES_APPROVAL"
BLOCKED_ORIGINAL = "BLOCKED_ORIGINAL"

ACCEPTED = "ACCEPTED"
ACCEPTED_WITH_CAVEAT = "ACCEPTED_WITH_CAVEAT"
NOT_ACCEPTED_FOR_EQUIVALENCE = "NOT_ACCEPTED_FOR_EQUIVALENCE"
SKIPPED = "SKIPPED"

SOURCE_TRACEABILITY_COLUMNS = [
    "course_module",
    "source_script",
    "original_key_calls",
    "original_logic_goal",
    "adapted_script_or_output",
    "adapted_key_calls",
    "logic_preserved",
    "deviation_from_original",
    "equivalence_level",
    "acceptance_status",
    "validation_evidence",
    "python_only_feasibility",
    "python_replacement_notes",
]


@dataclass(frozen=True)
class TraceRow:
    course_module: str
    source_script: str
    original_key_calls: str
    original_logic_goal: str
    adapted_script_or_output: str
    adapted_key_calls: str
    logic_preserved: str
    deviation_from_original: str
    equivalence_level: str
    acceptance_status: str
    evidence_files: tuple[str, ...]
    python_only_feasibility: str
    python_replacement_notes: str


def rel(*parts: str) -> str:
    return "/".join(parts)


TRACE_ROWS = (
    TraceRow(
        course_module="02_r_package_install",
        source_script=rel("2.R包安装", "1.基础R包安装.r")
        + ";"
        + rel("2.R包安装", "2.Github安装.r"),
        original_key_calls="install.packages(); BiocManager::install(); devtools::install_github(); local source package installation",
        original_logic_goal="Prepare the R package stack required by the course scripts.",
        adapted_script_or_output=".conda/seurat-core; .conda/scverse-course; sessionInfo logs",
        adapted_key_calls="conda environments; R package import checks; Python package import checks",
        logic_preserved="Records the dependency inventory and validates the runnable R/Python packages needed by the adapted workflow.",
        deviation_from_original="Does not blindly install every course package; unavailable Bioconductor/GitHub-heavy packages are blocked or replaced.",
        equivalence_level=TRACEABLE_REPLACEMENT,
        acceptance_status=ACCEPTED_WITH_CAVEAT,
        evidence_files=(
            "analysis/seurat_v5_course_run/logs/sessionInfo.txt",
            "analysis/seurat_v5_logic_run/tables/module_status.tsv",
        ),
        python_only_feasibility="PYTHON_EQUIVALENT_AVAILABLE",
        python_replacement_notes="A Python-only environment can be pinned with conda/pip, but exact R package execution would no longer be part of the run.",
    ),
    TraceRow(
        course_module="03_single_matrix_import",
        source_script=rel("3.单矩阵读取", "单矩阵读取.r"),
        original_key_calls="CreateSeuratObject(); PercentageFeatureSet(pattern='^MT-'); PercentageFeatureSet(pattern='^RP')",
        original_logic_goal="Read one expression matrix, create a Seurat object, and compute mitochondrial/ribosomal QC metrics.",
        adapted_script_or_output="scripts/run_seurat_v5_course_core.R; analysis/seurat_v5_course_run/objects/single_matrix_import.rds",
        adapted_key_calls="CreateSeuratObject(); PercentageFeatureSet()",
        logic_preserved="Uses the real course matrix and preserves Seurat object construction plus percent.mt/percent.rb metrics.",
        deviation_from_original="Explicit paths and saveRDS are used instead of interactive setwd/qs output.",
        equivalence_level=COURSE_CONSISTENT,
        acceptance_status=ACCEPTED,
        evidence_files=("analysis/seurat_v5_course_run/objects/single_matrix_import.rds",),
        python_only_feasibility="PYTHON_EQUIVALENT_AVAILABLE",
        python_replacement_notes="scanpy.read_text/read_mtx plus AnnData QC can reproduce the import logic, but it would not execute Seurat object creation.",
    ),
    TraceRow(
        course_module="04_multi_matrix_import",
        source_script=rel("4.多矩阵读取", "多矩阵读取.r"),
        original_key_calls="CreateSeuratObject(); merge(); PercentageFeatureSet()",
        original_logic_goal="Read multiple expression matrices and keep sample-specific metadata.",
        adapted_script_or_output="scripts/run_seurat_v5_course_core.R; multi_matrix_GSE118389.rds; multi_matrix_GSE158631.rds",
        adapted_key_calls="CreateSeuratObject(); PercentageFeatureSet()",
        logic_preserved="Imports the two real course matrix files and records per-object evidence.",
        deviation_from_original="The current validated run stores the imported objects separately; merge mechanics are validated in the multi-10X branch.",
        equivalence_level=COURSE_CONSISTENT,
        acceptance_status=ACCEPTED,
        evidence_files=(
            "analysis/seurat_v5_course_run/objects/multi_matrix_GSE118389.rds",
            "analysis/seurat_v5_course_run/objects/multi_matrix_GSE158631.rds",
        ),
        python_only_feasibility="PYTHON_EQUIVALENT_AVAILABLE",
        python_replacement_notes="AnnData concatenation can preserve the same sample labels, with different object class semantics.",
    ),
    TraceRow(
        course_module="05_standard_10x_import",
        source_script=rel("5.10X数据读取(标准10X)", "10X数据读取.r"),
        original_key_calls="Read10X(); CreateSeuratObject(); PercentageFeatureSet()",
        original_logic_goal="Read standard 10X folders from the course tar archive.",
        adapted_script_or_output="scripts/run_seurat_v5_course_core.R; standard10x_tar_import_list.rds",
        adapted_key_calls="Read10X(); CreateSeuratObject()",
        logic_preserved="Normalizes the real course standard 10X tar layout and reads it with Seurat.",
        deviation_from_original="Uses explicit extracted directories and saveRDS rather than interactive working directories.",
        equivalence_level=COURSE_CONSISTENT,
        acceptance_status=ACCEPTED,
        evidence_files=("analysis/seurat_v5_course_run/objects/standard10x_tar_import_list.rds",),
        python_only_feasibility="PYTHON_EQUIVALENT_AVAILABLE",
        python_replacement_notes="scanpy.read_10x_mtx can reproduce 10X import and QC metrics.",
    ),
    TraceRow(
        course_module="06_nonstandard_10x_import",
        source_script=rel("6.10X数据读取(非标准10X)", "非标准10X数据读取 .R"),
        original_key_calls="R.utils decompression; Read10X-compatible layout; CreateSeuratObject()",
        original_logic_goal="Normalize a non-standard 10X archive layout before object creation.",
        adapted_script_or_output="scripts/run_seurat_v5_course_core.R; nonstandard10x_first_import.rds",
        adapted_key_calls="file normalization; Read10X(); CreateSeuratObject()",
        logic_preserved="Uses the real nested course files and converts them to a readable 10X layout.",
        deviation_from_original="Conversion is scripted deterministically instead of relying on manual file operations.",
        equivalence_level=COURSE_CONSISTENT,
        acceptance_status=ACCEPTED,
        evidence_files=("analysis/seurat_v5_course_run/objects/nonstandard10x_first_import.rds",),
        python_only_feasibility="PYTHON_EQUIVALENT_AVAILABLE",
        python_replacement_notes="Python can rename/decompress files and read the resulting Matrix Market files into AnnData.",
    ),
    TraceRow(
        course_module="07_multiple_10x_merge",
        source_script=rel("7.多个10X数据读取", "多个10X数据读取.r"),
        original_key_calls="Read10X(); CreateSeuratObject(); merge()",
        original_logic_goal="Read multiple 10X samples and merge them while preserving sample identity.",
        adapted_script_or_output="scripts/run_seurat_v5_course_core.R; multi10x_merged_raw.rds",
        adapted_key_calls="Read10X(); CreateSeuratObject(); merge()",
        logic_preserved="Reads the real GSE184198 course samples and merges them with sample_id/condition metadata.",
        deviation_from_original="Explicit sample metadata and output paths replace interactive course paths.",
        equivalence_level=COURSE_CONSISTENT,
        acceptance_status=ACCEPTED,
        evidence_files=("analysis/seurat_v5_course_run/objects/multi10x_merged_raw.rds",),
        python_only_feasibility="PYTHON_EQUIVALENT_AVAILABLE",
        python_replacement_notes="scanpy.concat can reproduce the merge logic and sample metadata preservation.",
    ),
    TraceRow(
        course_module="08_h5_import",
        source_script=rel("8.H5数据读取", "H5读取.r"),
        original_key_calls="Read10X_h5(); CreateSeuratObject(); PercentageFeatureSet()",
        original_logic_goal="Read 10X H5 matrices and build Seurat objects.",
        adapted_script_or_output="scripts/run_seurat_v5_course_core.R; h5_first_import.rds",
        adapted_key_calls="Read10X_h5(); CreateSeuratObject()",
        logic_preserved="Extracts real course H5 files and validates at least one H5 import path.",
        deviation_from_original="Only the first H5 import product is retained as a smoke-tested object.",
        equivalence_level=COURSE_CONSISTENT,
        acceptance_status=ACCEPTED_WITH_CAVEAT,
        evidence_files=("analysis/seurat_v5_course_run/objects/h5_first_import.rds",),
        python_only_feasibility="PYTHON_EQUIVALENT_AVAILABLE",
        python_replacement_notes="scanpy.read_10x_h5 can replace the H5 import logic.",
    ),
    TraceRow(
        course_module="09_mixed_format_merge",
        source_script=rel("9.多个不同数据合并(思路)", "多个不同数据合并(思路).r"),
        original_key_calls="Seurat object import per source; metadata harmonization; merge()",
        original_logic_goal="Describe how to merge different single-cell input formats after converting each to a Seurat object.",
        adapted_script_or_output="scripts/run_seurat_v5_course_core.R; multi10x_merged_raw.rds",
        adapted_key_calls="format-specific import branches; merge(); sample metadata preservation",
        logic_preserved="The validated run demonstrates the same merge principle on real course 10X samples with harmonized metadata.",
        deviation_from_original="The original script is a conceptual branch; no additional mixed-format data product is generated beyond validated merge mechanics.",
        equivalence_level=TRACEABLE_REPLACEMENT,
        acceptance_status=ACCEPTED_WITH_CAVEAT,
        evidence_files=("analysis/seurat_v5_course_run/objects/multi10x_merged_raw.rds",),
        python_only_feasibility="PYTHON_EQUIVALENT_AVAILABLE",
        python_replacement_notes="Python can read each format into AnnData and concatenate after metadata harmonization.",
    ),
    TraceRow(
        course_module="10_qc",
        source_script=rel("10.数据质控", "数据质控.r"),
        original_key_calls="VlnPlot(); FeatureScatter(); subset(nCount_RNA,nFeature_RNA,percent.mt,percent.rb)",
        original_logic_goal="Visualize QC metrics and filter low-quality cells.",
        adapted_script_or_output="scripts/run_seurat_v5_course_core.R; qc_cell_counts.tsv; QC figures",
        adapted_key_calls="VlnPlot(); subset(); qc count table",
        logic_preserved="Runs pre/post QC plots and records cell counts on the real merged course object.",
        deviation_from_original="Uses deterministic thresholds selected for the validated course data rather than interactive inspection.",
        equivalence_level=COURSE_CONSISTENT,
        acceptance_status=ACCEPTED,
        evidence_files=(
            "analysis/seurat_v5_course_run/tables/qc_cell_counts.tsv",
            "analysis/seurat_v5_course_run/figures/01_QC_before_violin.png",
            "analysis/seurat_v5_course_run/figures/02_QC_after_violin.png",
        ),
        python_only_feasibility="PYTHON_EQUIVALENT_AVAILABLE",
        python_replacement_notes="scanpy.pp.calculate_qc_metrics plus seaborn/matplotlib can reproduce QC plots and filters.",
    ),
    TraceRow(
        course_module="11_normalization_harmony",
        source_script=rel("11.标准化(过滤游离RNA)", "标准化步骤.r"),
        original_key_calls="NormalizeData(); FindVariableFeatures(); ScaleData(); RunPCA(); RunTSNE(); RunUMAP(); decontX(); RunHarmony()",
        original_logic_goal="Normalize expression, reduce dimensions, optionally decontaminate counts, and integrate batches with Harmony.",
        adapted_script_or_output="scripts/run_seurat_v5_course_core.R; processed_multi10x_qc_harmony_cluster.rds",
        adapted_key_calls="NormalizeData(); FindVariableFeatures(); ScaleData(); RunPCA(); RunHarmony(); RunUMAP()",
        logic_preserved="Core Seurat normalization, variable features, PCA, Harmony integration, and UMAP are run on real course data.",
        deviation_from_original="decontX branch is not executed because celda/decontX dependencies were unavailable.",
        equivalence_level=METHOD_EQUIVALENT,
        acceptance_status=ACCEPTED_WITH_CAVEAT,
        evidence_files=(
            "analysis/seurat_v5_course_run/objects/processed_multi10x_qc_harmony_cluster.rds",
            "analysis/seurat_v5_course_run/figures/03_variable_features.png",
            "analysis/seurat_v5_course_run/figures/04_elbow_plot.png",
        ),
        python_only_feasibility="PYTHON_EQUIVALENT_AVAILABLE",
        python_replacement_notes="scanpy normalization/PCA/UMAP plus harmonypy can reproduce the analysis intent; decontX would need a Python ambient-RNA alternative.",
    ),
    TraceRow(
        course_module="12_doublet",
        source_script=rel("12-1.去除双细胞(DoubletFinder)", "DoubletFinder.r")
        + ";"
        + rel("12-2.去除双细胞(scDblFinder)", "scDblFinder.r"),
        original_key_calls="paramSweep(); find.pK(); modelHomotypic(); doubletFinder(); scDblFinder(samples=...)",
        original_logic_goal="Detect likely doublets per sample and subset singlets.",
        adapted_script_or_output="scripts/run_logic_scrublet.py; scrublet_doublet_scores.tsv; cell_counts_before_after_doublet.tsv",
        adapted_key_calls="scrublet.Scrublet(); per-sample QC top-5-percent fallback",
        logic_preserved="Produces per-sample doublet calls, removes predicted doublets, and carries singlets into reprocessing.",
        deviation_from_original="Neither DoubletFinder nor scDblFinder ran; Scrublet automatic threshold failed and a documented QC fallback was used.",
        equivalence_level=TRACEABLE_REPLACEMENT,
        acceptance_status=ACCEPTED_WITH_CAVEAT,
        evidence_files=(
            "analysis/seurat_v5_logic_run/tables/scrublet_doublet_scores.tsv",
            "analysis/seurat_v5_logic_run/tables/cell_counts_before_after_doublet.tsv",
            "analysis/seurat_v5_logic_run/figures/umap_doublet_class.png",
        ),
        python_only_feasibility="PYTHON_EQUIVALENT_AVAILABLE",
        python_replacement_notes="Scrublet is a standard Python doublet detector, but it is not numerically equivalent to DoubletFinder/scDblFinder.",
    ),
    TraceRow(
        course_module="12_3_post_doublet",
        source_script=rel("12-3.去除双细胞后常规标准化", "去除双细胞后常规标准化.r"),
        original_key_calls="NormalizeData(); FindVariableFeatures(); ScaleData(); RunHarmony(); RunUMAP(); RunTSNE()",
        original_logic_goal="Re-run the standard Seurat workflow after doublet removal.",
        adapted_script_or_output="scripts/run_logic_reprocess_singlets.R; annotated_consensus.rds",
        adapted_key_calls="NormalizeData(); FindVariableFeatures(); ScaleData(); RunPCA(); RunHarmony(); RunUMAP(); FindAllMarkers()",
        logic_preserved="Uses the singlet set and reprocesses normalization, PCA, Harmony, UMAP, clustering, markers, and annotation.",
        deviation_from_original="The singlet set comes from the replacement doublet workflow rather than DoubletFinder/scDblFinder.",
        equivalence_level=METHOD_EQUIVALENT,
        acceptance_status=ACCEPTED_WITH_CAVEAT,
        evidence_files=(
            "analysis/seurat_v5_logic_run/objects/rds/annotated_consensus.rds",
            "analysis/seurat_v5_logic_run/tables/singlet_reprocess_summary.tsv",
        ),
        python_only_feasibility="PYTHON_EQUIVALENT_AVAILABLE",
        python_replacement_notes="A Scanpy-only reprocessing path is feasible, but current accepted object is a Seurat RDS.",
    ),
    TraceRow(
        course_module="13_clustering",
        source_script=rel("13.分群聚类", "分群聚类.r"),
        original_key_calls="FindNeighbors(reduction='harmony'); FindClusters(); clustree; RunUMAP(); RunTSNE()",
        original_logic_goal="Cluster cells after dimensional reduction and visualize embeddings.",
        adapted_script_or_output="scripts/run_seurat_v5_course_core.R; processed object and UMAP figures",
        adapted_key_calls="FindNeighbors(); FindClusters(); RunUMAP()",
        logic_preserved="Clusters the real course object on Harmony embeddings and writes UMAP figures.",
        deviation_from_original="Resolution sweep/clustree is represented by validated selected clustering rather than all interactive plots.",
        equivalence_level=METHOD_EQUIVALENT,
        acceptance_status=ACCEPTED_WITH_CAVEAT,
        evidence_files=(
            "analysis/seurat_v5_course_run/figures/05_umap_clusters.png",
            "analysis/seurat_v5_course_run/figures/06_umap_sample_id.png",
            "analysis/seurat_v5_course_run/figures/07_umap_condition.png",
        ),
        python_only_feasibility="PYTHON_EQUIVALENT_AVAILABLE",
        python_replacement_notes="scanpy.pp.neighbors and scanpy.tl.leiden can reproduce the clustering intent but not Seurat cluster IDs exactly.",
    ),
    TraceRow(
        course_module="14_manual_annotation",
        source_script=rel("14.细胞手动注释", "细胞手动注释.r"),
        original_key_calls="DotPlot(); FeaturePlot(); DimPlot(); subcluster immune cells; assign cellType metadata",
        original_logic_goal="Assign biologically supported cell-type labels from marker expression and embeddings.",
        adapted_script_or_output="scripts/run_logic_reprocess_singlets.R; consensus_annotation.tsv; marker_rule_scores.tsv",
        adapted_key_calls="FindAllMarkers(); curated marker scoring; UMAP marker-rule plots",
        logic_preserved="Uses marker evidence to assign broad consensus labels and writes marker-score evidence per cell.",
        deviation_from_original="Labels are broad marker-rule labels, not manual expert-curated labels from the exact course marker panels.",
        equivalence_level=TRACEABLE_REPLACEMENT,
        acceptance_status=ACCEPTED_WITH_CAVEAT,
        evidence_files=(
            "analysis/seurat_v5_logic_run/tables/consensus_annotation.tsv",
            "analysis/seurat_v5_logic_run/tables/marker_rule_scores.tsv",
            "analysis/seurat_v5_logic_run/figures/marker_dotplot_consensus.png",
        ),
        python_only_feasibility="PYTHON_EQUIVALENT_AVAILABLE",
        python_replacement_notes="Python marker scoring and plotting are feasible; expert label validation remains a biological review step.",
    ),
    TraceRow(
        course_module="27_subset_extraction",
        source_script=rel("14.细胞手动注释", "细胞手动注释.r"),
        original_key_calls="subset(); immune compartment extraction; re-normalize/recluster selected populations",
        original_logic_goal="Extract a biologically defined cell population for focused downstream analysis.",
        adapted_script_or_output="scripts/run_logic_subset.R; subset_T_cells.rds; subset_T_cells_barcodes.tsv",
        adapted_key_calls="subset(obj, cells=...); saveRDS(); barcode membership and summary tables",
        logic_preserved="Creates a reusable Seurat subset object from annotated cells and records the exact subset rule.",
        deviation_from_original="The default subset is parameterized T_cells rather than the exact immune/tumor subsetting choices made interactively in the course script.",
        equivalence_level=COURSE_CONSISTENT,
        acceptance_status=ACCEPTED_WITH_CAVEAT,
        evidence_files=(
            "analysis/seurat_v5_logic_run/objects/rds/subset_T_cells.rds",
            "analysis/seurat_v5_logic_run/tables/subset_T_cells_barcodes.tsv",
            "analysis/seurat_v5_logic_run/tables/subset_T_cells_summary.tsv",
        ),
        python_only_feasibility="PYTHON_EQUIVALENT_AVAILABLE",
        python_replacement_notes="AnnData subsetting can reproduce the data operation, but Seurat RDS subset output requires R.",
    ),
    TraceRow(
        course_module="29_gene_expression_plotting",
        source_script=rel("14.细胞手动注释", "细胞手动注释.r"),
        original_key_calls="FeaturePlot(); VlnPlot(); DotPlot(); DimPlot() marker visualization",
        original_logic_goal="Display marker or candidate gene expression on embeddings and grouped summaries.",
        adapted_script_or_output="scripts/run_logic_gene_expression_plotting.R; gene_expression_* tables and figures",
        adapted_key_calls="FeaturePlot(); VlnPlot(); DotPlot(); FetchData(layer=...)",
        logic_preserved="Uses the final annotated Seurat object to plot requested genes on UMAP and summarize expression by cell type.",
        deviation_from_original="The gene panel is a parameterized default marker panel, not the exact interactive marker lists from the course.",
        equivalence_level=COURSE_CONSISTENT,
        acceptance_status=ACCEPTED_WITH_CAVEAT,
        evidence_files=(
            "analysis/seurat_v5_logic_run/tables/gene_expression_requested_genes.tsv",
            "analysis/seurat_v5_logic_run/tables/gene_expression_summary.tsv",
            "analysis/seurat_v5_logic_run/figures/gene_expression_featureplot.png",
            "analysis/seurat_v5_logic_run/figures/gene_expression_dotplot.png",
        ),
        python_only_feasibility="PYTHON_EQUIVALENT_AVAILABLE",
        python_replacement_notes="Scanpy/matplotlib can reproduce expression plots, but the current accepted implementation intentionally uses Seurat plotting.",
    ),
    TraceRow(
        course_module="15_16_auto_annotation",
        source_script=rel("15-1.细胞自动注释(singleR)", "single.r")
        + ";"
        + rel("15-2.细胞自动注释(SCINA)", "SCINA.r")
        + ";"
        + rel("16-1.自有模型细胞自动注释(TransferData)", "TransferData.R")
        + ";"
        + rel("16-2.自有模型细胞自动注释(scPred)", "scPred.r"),
        original_key_calls="SingleR(); SCINA(); FindTransferAnchors(); TransferData(); getFeatureSpace(); trainModel(); scPredict()",
        original_logic_goal="Annotate cells using external references, marker signatures, or trained classifiers.",
        adapted_script_or_output="scripts/run_logic_celltypist_bounded.py; celltypist_bounded_annotation.tsv; consensus_annotation.tsv",
        adapted_key_calls="local CellTypist model resolution; marker-rule consensus fallback",
        logic_preserved="Keeps automated/reference annotation as secondary evidence and validates labels against marker-rule consensus.",
        deviation_from_original="SingleR, SCINA, TransferData, and scPred are not executed; CellTypist is skipped unless a local .pkl model exists.",
        equivalence_level=SKIPPED_DEPENDENCY,
        acceptance_status=SKIPPED,
        evidence_files=(
            "analysis/seurat_v5_logic_run/tables/celltypist_bounded_annotation.tsv",
            "analysis/seurat_v5_logic_run/tables/consensus_annotation.tsv",
        ),
        python_only_feasibility="PYTHON_EQUIVALENT_AVAILABLE",
        python_replacement_notes="CellTypist and scikit-learn classifiers can provide Python annotation, but exact SingleR/SCINA/TransferData/scPred equivalence is not claimed.",
    ),
    TraceRow(
        course_module="15_llm_annotation",
        source_script=rel("15-3.细胞自动注释(调用Kimi)", "调用Kimi进行自动注释.r")
        + ";"
        + rel("15-4.细胞自动注释(调用DeepSeek)", "调用DeepSeek进行自动注释.r"),
        original_key_calls="FindAllMarkers(); xxdAIcelltype/openai client; Kimi/DeepSeek-compatible API",
        original_logic_goal="Send marker summaries to an LLM API for suggested cell-type labels.",
        adapted_script_or_output="No external LLM call; local marker-rule consensus only",
        adapted_key_calls="No API call",
        logic_preserved="Preserves the privacy guardrail by using local marker summaries and not uploading data.",
        deviation_from_original="LLM-assisted labels are not produced because API approval/credentials were not provided.",
        equivalence_level=SKIPPED_REQUIRES_APPROVAL,
        acceptance_status=SKIPPED,
        evidence_files=("analysis/seurat_v5_logic_run/tables/consensus_annotation.tsv",),
        python_only_feasibility="PYTHON_EQUIVALENT_AVAILABLE",
        python_replacement_notes="Python API clients could reproduce the LLM call only after explicit approval and credentials.",
    ),
    TraceRow(
        course_module="17_marker_algorithms",
        source_script=rel("17.多种算法计算marker基因", "多种算法计算marker基因.R"),
        original_key_calls="FindAllMarkers(); presto::wilcoxauc(); COSG::cosg(); starTracer::searchMarker()",
        original_logic_goal="Compute marker genes with multiple algorithms and compare evidence.",
        adapted_script_or_output="scripts/run_logic_reprocess_singlets.R; cluster_markers_seurat.tsv",
        adapted_key_calls="FindAllMarkers()",
        logic_preserved="Runs Seurat marker detection on the final singlet object and writes top marker evidence.",
        deviation_from_original="presto, COSG, and starTracer branches are not executed.",
        equivalence_level=METHOD_EQUIVALENT,
        acceptance_status=ACCEPTED_WITH_CAVEAT,
        evidence_files=(
            "analysis/seurat_v5_logic_run/tables/cluster_markers_seurat.tsv",
            "analysis/seurat_v5_logic_run/tables/cluster_markers_top5_seurat.tsv",
        ),
        python_only_feasibility="PYTHON_EQUIVALENT_AVAILABLE",
        python_replacement_notes="scanpy.tl.rank_genes_groups can replace Seurat markers; exact FindAllMarkers statistics differ.",
    ),
    TraceRow(
        course_module="28_findmarkers_pairwise",
        source_script=rel("17.多种算法计算marker基因", "多种算法计算marker基因.R"),
        original_key_calls="FindAllMarkers(); marker comparisons using Seurat differential expression statistics",
        original_logic_goal="Find genes that distinguish cell groups and support annotation or comparison.",
        adapted_script_or_output="scripts/run_logic_findmarkers.R; findmarkers_pairwise.tsv; findmarkers_pairwise_volcano.png",
        adapted_key_calls="FindMarkers(ident.1=..., ident.2=..., group.by-derived identities)",
        logic_preserved="Uses Seurat differential expression on the final annotated object for a parameterized two-group comparison.",
        deviation_from_original="Runs one explicit pairwise comparison by default rather than all marker algorithms in the course marker script.",
        equivalence_level=COURSE_CONSISTENT,
        acceptance_status=ACCEPTED,
        evidence_files=(
            "analysis/seurat_v5_logic_run/tables/findmarkers_pairwise.tsv",
            "analysis/seurat_v5_logic_run/tables/findmarkers_pairwise_params.tsv",
            "analysis/seurat_v5_logic_run/figures/findmarkers_pairwise_volcano.png",
        ),
        python_only_feasibility="PYTHON_EQUIVALENT_AVAILABLE",
        python_replacement_notes="Scanpy Wilcoxon can replace the comparison logic, but the accepted module uses Seurat FindMarkers for course consistency.",
    ),
    TraceRow(
        course_module="18_go_kegg",
        source_script=rel("18.GO_KEGG", "GO_KEGG.r"),
        original_key_calls="clusterProfiler::enrichGO(); enrichKEGG(); org.Hs.eg.db/bitr()",
        original_logic_goal="Interpret marker gene lists through GO and KEGG enrichment.",
        adapted_script_or_output="scripts/run_logic_enrichment.py; enrichment_gprofiler.tsv",
        adapted_key_calls="gprofiler-official with request timeout when available; local marker-overlap fallback",
        logic_preserved="Uses marker gene lists for functional interpretation and writes an enrichment-like table.",
        deviation_from_original="clusterProfiler GO/KEGG did not run; online gProfiler is attempted with a bounded timeout and local overlap fallback is used when unavailable.",
        equivalence_level=TRACEABLE_REPLACEMENT,
        acceptance_status=ACCEPTED_WITH_CAVEAT,
        evidence_files=(
            "analysis/seurat_v5_logic_run/tables/enrichment_gprofiler.tsv",
            "analysis/seurat_v5_logic_run/tables/enrichment_method_status.tsv",
            "analysis/seurat_v5_logic_run/figures/enrichment_dotplot.png",
        ),
        python_only_feasibility="PYTHON_EQUIVALENT_AVAILABLE",
        python_replacement_notes="gprofiler-official or goatools can implement Python enrichment; online availability and database versions must be recorded.",
    ),
    TraceRow(
        course_module="19_pseudotime",
        source_script=rel("19-1.拟时序分析Monocle2", "拟时序分析Monocle2.r")
        + ";"
        + rel("19-2.拟时序分析Monocle3", "Monocle3.r"),
        original_key_calls="FindAllMarkers(); newCellDataSet/as.cell_data_set(); learn_graph(); order_cells()",
        original_logic_goal="Order a biologically plausible cell subset along a trajectory/pseudotime axis.",
        adapted_script_or_output="scripts/run_logic_trajectory.py; trajectory_pseudotime.tsv",
        adapted_key_calls="scanpy neighbors; diffusion map/DPT-style pseudotime; PAGA when possible",
        logic_preserved="Selects a large coherent compartment, computes graph-based pseudotime, and overlays pseudotime on UMAP.",
        deviation_from_original="Monocle2/3 graph learning and root selection are not numerically reproduced.",
        equivalence_level=TRACEABLE_REPLACEMENT,
        acceptance_status=ACCEPTED_WITH_CAVEAT,
        evidence_files=(
            "analysis/seurat_v5_logic_run/tables/trajectory_pseudotime.tsv",
            "analysis/seurat_v5_logic_run/figures/trajectory_pseudotime_umap.png",
        ),
        python_only_feasibility="PYTHON_EQUIVALENT_AVAILABLE",
        python_replacement_notes="Scanpy, Palantir, or scVelo can provide Python trajectory analyses, but method identity changes.",
    ),
    TraceRow(
        course_module="20_cellchat",
        source_script=rel("20.CellChat", "CellChat.r"),
        original_key_calls="createCellChat(); CellChatDB.human; subsetDB(); computeCommunProb(); netAnalysis_computeCentrality()",
        original_logic_goal="Infer cell-cell communication hypotheses from cell-type expression and ligand-receptor databases.",
        adapted_script_or_output="scripts/run_logic_communication.py; communication_lr_scores.tsv",
        adapted_key_calls="curated ligand-receptor mean-expression scoring",
        logic_preserved="Uses consensus cell types and ligand/receptor expression to rank sender-receiver hypotheses.",
        deviation_from_original="No CellChat probability model, pathway aggregation, centrality, or CellChatDB full database is reproduced.",
        equivalence_level=PROXY_ONLY,
        acceptance_status=NOT_ACCEPTED_FOR_EQUIVALENCE,
        evidence_files=(
            "analysis/seurat_v5_logic_run/tables/communication_lr_scores.tsv",
            "analysis/seurat_v5_logic_run/figures/communication_lr_heatmap.png",
        ),
        python_only_feasibility="PYTHON_PROXY_ONLY",
        python_replacement_notes="Python LIANA/CellPhoneDB-style tools could improve method fidelity; current implementation is only a local CellChat-like proxy.",
    ),
    TraceRow(
        course_module="21_22_cnv",
        source_script=rel("21.copykat", "copykat.r")
        + ";"
        + rel("22-1.infercnv", "infercnv.r")
        + ";"
        + rel("22-2.infercnv(引入正常对照)", "infercnv(引入正常对照).r")
        + ";"
        + rel("22-3.infercnv(计算打分)", "infercnv(计算打分).r"),
        original_key_calls="copykat(rawmat=...); CreateInfercnvObject(); infercnv::run(); inferCNV score parsing",
        original_logic_goal="Infer malignant/diploid or CNV-like expression programs using raw counts and reference groups.",
        adapted_script_or_output="scripts/run_logic_cnv_proxy.py; cnv_proxy_scores.tsv",
        adapted_key_calls="high-variance expression deviation score for epithelial/tumor-like cells versus reference cells",
        logic_preserved="Focuses on tumor-like epithelial cells and compares them to non-epithelial reference expression.",
        deviation_from_original="No copykat model, gene-order inferCNV heatmap, or validated normal-control reference run is reproduced.",
        equivalence_level=PROXY_ONLY,
        acceptance_status=NOT_ACCEPTED_FOR_EQUIVALENCE,
        evidence_files=(
            "analysis/seurat_v5_logic_run/tables/cnv_proxy_scores.tsv",
            "analysis/seurat_v5_logic_run/figures/cnv_proxy_umap.png",
        ),
        python_only_feasibility="PYTHON_PROXY_ONLY",
        python_replacement_notes="infercnvpy could provide a closer Python replacement if gene-order annotations and reference groups are fixed.",
    ),
    TraceRow(
        course_module="23_hdwgcna",
        source_script=rel("23.hdWGCNA", "hdWGCNA.r"),
        original_key_calls="SetupForWGCNA(); MetacellsByGroups(); NormalizeMetacells(); TestSoftPowers(); ConstructNetwork(); ModuleEigengenes(); ModuleFeaturePlot()",
        original_logic_goal="Build metacell-based co-expression modules and score module programs by cell type.",
        adapted_script_or_output="scripts/run_logic_coexpression.py; coexpression_gene_modules.tsv; coexpression_module_scores.tsv",
        adapted_key_calls="high-variance genes; metacell profiles; gene-gene correlation; hierarchical modules; module scores",
        logic_preserved="Preserves metacell construction, co-expression module detection, module membership, and module-celltype scoring.",
        deviation_from_original="No WGCNA soft-threshold selection, topological overlap network, or hdWGCNA hub-gene workflow is reproduced.",
        equivalence_level=TRACEABLE_REPLACEMENT,
        acceptance_status=ACCEPTED_WITH_CAVEAT,
        evidence_files=(
            "analysis/seurat_v5_logic_run/tables/coexpression_gene_modules.tsv",
            "analysis/seurat_v5_logic_run/tables/coexpression_module_scores.tsv",
            "analysis/seurat_v5_logic_run/tables/coexpression_module_celltype_summary.tsv",
            "analysis/seurat_v5_logic_run/figures/coexpression_module_celltype_heatmap.png",
        ),
        python_only_feasibility="PYTHON_EQUIVALENT_AVAILABLE",
        python_replacement_notes="PyWGCNA or a custom metacell/correlation pipeline can reproduce the analysis intent; hdWGCNA-specific outputs are not exact.",
    ),
    TraceRow(
        course_module="24_1_cibersort",
        source_script=rel("24-1反卷积(CIBERSORT)", "反卷积(CIBERSORT).R")
        + ";"
        + rel("24-1反卷积(CIBERSORT)", "xxdimmune.R"),
        original_key_calls="CIBERSORT(signature, mixture, perm=2, QN=F); e1071/limma/preprocessCore helper",
        original_logic_goal="Estimate bulk immune fractions from a signature matrix and bulk mixture expression.",
        adapted_script_or_output="scripts/run_seurat_v5_course_cibersort.R; CIBERSORT-Results.tsv",
        adapted_key_calls="source xxdimmune.R; CIBERSORT('scimmune.txt','uniq.symbol.txt',perm=2,QN=F)",
        logic_preserved="Uses the course helper, course TCGA expression, course LM genes, and processed Seurat-derived signature.",
        deviation_from_original="Runs a bounded 30-sample subset for feasible validation.",
        equivalence_level=COURSE_CONSISTENT,
        acceptance_status=ACCEPTED_WITH_CAVEAT,
        evidence_files=(
            "analysis/seurat_v5_course_run/cibersort/CIBERSORT-Results.tsv",
            "analysis/seurat_v5_course_run/cibersort/cibersort_fraction_barplot.png",
        ),
        python_only_feasibility="PYTHON_PROXY_ONLY",
        python_replacement_notes="A Python NNLS/SVR deconvolution could mimic the goal, but exact bundled CIBERSORT helper execution is R code.",
    ),
    TraceRow(
        course_module="24_2_music",
        source_script=rel("24-2反卷积(MusiC)", "反卷积(MusiC).r"),
        original_key_calls="SingleCellExperiment(); MuSiC::music_prop()",
        original_logic_goal="Use single-cell reference labels to deconvolve bulk expression with MuSiC.",
        adapted_script_or_output="Not executed; CIBERSORT branch is the validated deconvolution result",
        adapted_key_calls="No MuSiC replacement executed",
        logic_preserved="The broader deconvolution theme is represented by CIBERSORT, not MuSiC.",
        deviation_from_original="MuSiC package and required annotated single-cell reference were unavailable.",
        equivalence_level=BLOCKED_ORIGINAL,
        acceptance_status=SKIPPED,
        evidence_files=("analysis/seurat_v5_course_run/cibersort/CIBERSORT-Results.tsv",),
        python_only_feasibility="PYTHON_PROXY_ONLY",
        python_replacement_notes="Python deconvolution alternatives exist, but they would not be MuSiC-equivalent without reimplementing its model.",
    ),
)


def read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t")


def resolve_evidence_path(item: str, logic_dir: Path) -> Path:
    default_prefix = Path("analysis/seurat_v5_logic_run")
    item_path = Path(item)
    try:
        relative = item_path.relative_to(default_prefix)
        return logic_dir / relative
    except ValueError:
        return item_path


def evidence_status(path: Path) -> str:
    if not path.exists():
        return "MISSING"
    if not path.is_file():
        return "NOT_FILE"
    return f"OK {path.stat().st_size} bytes"


def file_evidence(paths: tuple[str, ...], logic_dir: Path) -> str:
    parts = []
    for item in paths:
        path = resolve_evidence_path(item, logic_dir)
        parts.append(f"{path} ({evidence_status(path)})")
    return "; ".join(parts)


def invalid_source_paths(course_root: Path, source_script: str) -> list[str]:
    invalid = []
    course_root_resolved = course_root.resolve()
    for rel_path in source_script.split(";"):
        rel_path = rel_path.strip()
        if not rel_path:
            continue
        candidate = (course_root / rel_path).resolve()
        try:
            candidate.relative_to(course_root_resolved)
        except ValueError:
            invalid.append(f"{rel_path} escapes course root")
            continue
        if not candidate.exists():
            invalid.append(f"{rel_path} missing")
        elif not candidate.is_file():
            invalid.append(f"{rel_path} is not a file")
    return invalid


def build_traceability(course_root: Path, logic_dir: Path) -> pd.DataFrame:
    rows = []
    for row in TRACE_ROWS:
        rows.append(
            {
                "course_module": row.course_module,
                "source_script": row.source_script,
                "original_key_calls": row.original_key_calls,
                "original_logic_goal": row.original_logic_goal,
                "adapted_script_or_output": row.adapted_script_or_output,
                "adapted_key_calls": row.adapted_key_calls,
                "logic_preserved": row.logic_preserved,
                "deviation_from_original": row.deviation_from_original,
                "equivalence_level": row.equivalence_level,
                "acceptance_status": row.acceptance_status,
                "validation_evidence": file_evidence(row.evidence_files, logic_dir),
                "python_only_feasibility": row.python_only_feasibility,
                "python_replacement_notes": row.python_replacement_notes,
            }
        )
    return pd.DataFrame(rows, columns=SOURCE_TRACEABILITY_COLUMNS)


def validate_traceability(
    table: pd.DataFrame,
    course_root: Path,
    logic_dir: Path | None = None,
    *,
    require_evidence: bool = True,
) -> list[str]:
    _ = logic_dir
    errors: list[str] = []
    if list(table.columns) != SOURCE_TRACEABILITY_COLUMNS:
        errors.append("columns do not match SOURCE_TRACEABILITY_COLUMNS")

    required_columns = SOURCE_TRACEABILITY_COLUMNS
    for column in required_columns:
        if column not in table.columns:
            errors.append(f"missing required column: {column}")
            continue
        missing = table[column].fillna("").astype(str).str.strip().eq("")
        if missing.any():
            modules = ", ".join(table.loc[missing, "course_module"].astype(str).tolist())
            errors.append(f"{column} is blank for: {modules}")

    if "course_module" in table.columns and table["course_module"].duplicated().any():
        duplicated = table.loc[table["course_module"].duplicated(), "course_module"].astype(str).tolist()
        errors.append("duplicate course_module values: " + ", ".join(duplicated))

    allowed_levels = {
        COURSE_CONSISTENT,
        METHOD_EQUIVALENT,
        TRACEABLE_REPLACEMENT,
        PROXY_ONLY,
        SKIPPED_DEPENDENCY,
        SKIPPED_REQUIRES_APPROVAL,
        BLOCKED_ORIGINAL,
    }
    bad_levels = sorted(set(table.get("equivalence_level", pd.Series(dtype=str)).dropna()) - allowed_levels)
    if bad_levels:
        errors.append("unknown equivalence_level values: " + ", ".join(bad_levels))

    if {"equivalence_level", "acceptance_status"}.issubset(table.columns):
        bad_proxy = table[
            table["equivalence_level"].eq(PROXY_ONLY)
            & table["acceptance_status"].eq(ACCEPTED)
        ]
        if not bad_proxy.empty:
            errors.append(
                "PROXY_ONLY rows cannot be acceptance_status=ACCEPTED: "
                + ", ".join(bad_proxy["course_module"].astype(str).tolist())
            )

    if {"source_script", "course_module"}.issubset(table.columns):
        for _, row in table.iterrows():
            bad_sources = invalid_source_paths(course_root, str(row["source_script"]))
            if bad_sources:
                errors.append(
                    f"source_script is not a file for {row['course_module']}: "
                    + "; ".join(bad_sources)
                )

    if require_evidence and {"validation_evidence", "course_module", "acceptance_status"}.issubset(table.columns):
        for _, row in table.iterrows():
            evidence = str(row["validation_evidence"])
            if "MISSING" in evidence or "NOT_FILE" in evidence:
                errors.append(f"evidence missing for {row['course_module']}: {evidence}")

    return errors


def markdown_table(df: pd.DataFrame) -> str:
    if df.empty:
        return "_No rows._"
    lines = [
        "| " + " | ".join(df.columns) + " |",
        "| " + " | ".join(["---"] * len(df.columns)) + " |",
    ]
    for _, row in df.iterrows():
        values = [str(row[column]).replace("|", "/") for column in df.columns]
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def write_report(table: pd.DataFrame, report_path: Path) -> None:
    counts = (
        table.groupby(["equivalence_level", "acceptance_status"], dropna=False)
        .size()
        .reset_index(name="modules")
        .sort_values(["equivalence_level", "acceptance_status"])
    )
    proxy = table[table["equivalence_level"].eq(PROXY_ONLY)][
        ["course_module", "original_key_calls", "adapted_key_calls", "deviation_from_original", "acceptance_status"]
    ]
    python_summary = (
        table.groupby("python_only_feasibility", dropna=False)
        .size()
        .reset_index(name="modules")
        .sort_values("python_only_feasibility")
    )

    report = f"""# Seurat V5 Source Traceability Report

## Purpose

This report audits whether the current runnable workflow is traceable to the original Seurat V5 course scripts. It does not claim that replacement modules are numerically identical to unavailable original packages. Proxy modules are explicitly labeled and cannot be accepted as method-equivalent outputs.

## Equivalence Summary

{markdown_table(counts)}

## Python-Only Feasibility

{markdown_table(python_summary)}

## Proxy Modules That Are Not Method-Equivalent

{markdown_table(proxy)}

## Full Traceability Matrix

See `tables/source_traceability_matrix.tsv`.
"""
    report_path.write_text(report)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--course-root", required=True, type=Path)
    parser.add_argument("--logic-dir", required=True, type=Path)
    args = parser.parse_args()

    out_dir = args.logic_dir
    tables_dir = out_dir / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)

    table = build_traceability(args.course_root, args.logic_dir)
    errors = validate_traceability(table, args.course_root, args.logic_dir)
    if errors:
        for error in errors:
            print(f"ERROR: {error}")
        raise SystemExit(1)

    matrix_path = tables_dir / "source_traceability_matrix.tsv"
    report_path = out_dir / "source_traceability_report.md"
    table.to_csv(matrix_path, sep="\t", index=False)
    write_report(table, report_path)
    print(f"Wrote {matrix_path}")
    print(f"Wrote {report_path}")


if __name__ == "__main__":
    main()
