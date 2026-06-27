#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(harmony)
  library(ggplot2)
  library(dplyr)
})

parse_args <- function(args) {
  values <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--") || i == length(args)) {
      stop("Expected --key value arguments, got: ", paste(args, collapse = " "))
    }
    values[[sub("^--", "", key)]] <- args[[i + 1]]
    i <- i + 2
  }
  values
}

append_status <- function(out_dir, module, status, detail = "") {
  status_file <- file.path(out_dir, "tables", "module_status.tsv")
  if (!file.exists(status_file)) {
    writeLines("module\tstatus\tdetail", status_file)
  }
  line <- paste(module, status, gsub("[\t\r\n]+", " ", detail), sep = "\t")
  write(line, status_file, append = TRUE)
}

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

score_marker_sets <- function(obj, marker_sets) {
  expr <- GetAssayData(obj, assay = "RNA", layer = "data")
  cluster_ids <- levels(Idents(obj))
  rows <- list()

  for (cluster_id in cluster_ids) {
    cells <- WhichCells(obj, idents = cluster_id)
    for (label in names(marker_sets)) {
      present <- intersect(marker_sets[[label]], rownames(expr))
      score <- NA_real_
      if (length(present) > 0 && length(cells) > 0) {
        score <- mean(Matrix::colMeans(expr[present, cells, drop = FALSE]))
      }
      rows[[length(rows) + 1]] <- data.frame(
        cluster = cluster_id,
        candidate_celltype = label,
        marker_score = score,
        present_markers = paste(present, collapse = ","),
        stringsAsFactors = FALSE
      )
    }
  }

  scores <- dplyr::bind_rows(rows)
  best <- scores %>%
    group_by(cluster) %>%
    arrange(desc(marker_score), .by_group = TRUE) %>%
    mutate(rank = row_number()) %>%
    ungroup()

  top_specific <- best %>%
    filter(candidate_celltype != "Immune") %>%
    group_by(cluster) %>%
    arrange(desc(marker_score), .by_group = TRUE) %>%
    slice_head(n = 1) %>%
    ungroup()

  top_broad <- best %>%
    filter(candidate_celltype == "Immune") %>%
    select(cluster, broad_immune_score = marker_score, broad_immune_markers = present_markers)

  top <- top_specific %>%
    transmute(
      cluster,
      marker_rule_celltype = ifelse(is.na(marker_score) | marker_score < 0.02, "Unknown", candidate_celltype),
      top_marker_score = marker_score,
      top_present_markers = present_markers
    ) %>%
    left_join(top_broad, by = "cluster")

  second <- best %>%
    filter(rank == 2) %>%
    transmute(cluster, second_celltype = candidate_celltype, second_marker_score = marker_score)

  list(
    all_scores = scores,
    cluster_labels = left_join(top, second, by = "cluster")
  )
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
input_rds <- args[["input"]]
doublet_path <- args[["doublets"]]
out_dir <- args[["out"]]

if (is.null(input_rds) || is.null(doublet_path) || is.null(out_dir)) {
  stop("Usage: run_logic_reprocess_singlets.R --input <processed.rds> --doublets <doublet.tsv> --out <output_root>")
}
if (!file.exists(input_rds)) stop("Input RDS does not exist: ", input_rds)
if (!file.exists(doublet_path)) stop("Doublet TSV does not exist: ", doublet_path)

dir.create(file.path(out_dir, "objects", "rds"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "logs"), recursive = TRUE, showWarnings = FALSE)

set.seed(20260626)
append_status(out_dir, "post_doublet_reprocess", "START", input_rds)

obj <- readRDS(input_rds)
DefaultAssay(obj) <- "RNA"
if (inherits(obj[["RNA"]], "Assay5")) {
  obj <- JoinLayers(obj, assay = "RNA")
}

doublets <- read.delim(doublet_path, check.names = FALSE, stringsAsFactors = FALSE)
required_cols <- c("barcode", "doublet_score", "predicted_doublet", "method")
missing_cols <- setdiff(required_cols, colnames(doublets))
if (length(missing_cols) > 0) {
  stop("Doublet TSV missing required columns: ", paste(missing_cols, collapse = ", "))
}

match_idx <- match(colnames(obj), doublets$barcode)
if (anyNA(match_idx)) {
  stop("Doublet TSV missing ", sum(is.na(match_idx)), " cells from Seurat object")
}

predicted <- toupper(as.character(doublets$predicted_doublet[match_idx])) %in% c("TRUE", "T", "1", "YES")
obj$logic_doublet_score <- as.numeric(doublets$doublet_score[match_idx])
obj$logic_predicted_doublet <- predicted
obj$logic_doublet_method <- doublets$method[match_idx]

singlet_cells <- colnames(obj)[!predicted]
singlet <- subset(obj, cells = singlet_cells)
singlet@reductions <- list()
singlet@graphs <- list()
DefaultAssay(singlet) <- "RNA"
if (inherits(singlet[["RNA"]], "Assay5")) {
  singlet <- JoinLayers(singlet, assay = "RNA")
}

singlet <- NormalizeData(singlet, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
singlet <- FindVariableFeatures(singlet, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
singlet <- ScaleData(singlet, verbose = FALSE)
singlet <- RunPCA(singlet, features = VariableFeatures(singlet), npcs = 30, verbose = FALSE)
singlet <- harmony::RunHarmony(singlet, group.by.vars = "sample_id", reduction.use = "pca", dims.use = 1:20, verbose = FALSE)
singlet <- FindNeighbors(singlet, reduction = "harmony", dims = 1:20, verbose = FALSE)
singlet <- FindClusters(singlet, resolution = 0.5, algorithm = 1, verbose = FALSE)
singlet <- RunUMAP(singlet, reduction = "harmony", dims = 1:20, verbose = FALSE)
singlet <- JoinLayers(singlet, assay = "RNA")

saveRDS(singlet, file.path(out_dir, "objects", "rds", "singlet_processed_clustered.rds"))

Idents(singlet) <- "seurat_clusters"
markers <- FindAllMarkers(
  singlet,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  verbose = FALSE
)
write_tsv(markers, file.path(out_dir, "tables", "cluster_markers_seurat.tsv"))

top5 <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5, with_ties = FALSE) %>%
  ungroup()
write_tsv(top5, file.path(out_dir, "tables", "cluster_markers_top5_seurat.tsv"))

marker_sets <- list(
  Immune = c("PTPRC"),
  T_cells = c("CD3D", "CD3E", "CD2", "CD4", "CD8A", "NKG7", "GNLY"),
  B_cells = c("MS4A1", "CD79A", "CD19"),
  Plasma = c("MZB1", "SDC1", "IGHG1"),
  Myeloid_monocyte = c("LYZ", "CD14", "FCGR3A", "CD68", "CD163"),
  Epithelial_tumor_like = c("EPCAM", "KRT8", "KRT18", "KRT19"),
  Endothelial = c("PECAM1", "VWF", "KDR"),
  Fibroblast_stromal = c("COL1A1", "COL1A2", "DCN", "LUM", "MME"),
  Mast = c("TPSAB1", "TPSB2", "CPA3", "KIT")
)

marker_scores <- score_marker_sets(singlet, marker_sets)
write_tsv(marker_scores$all_scores, file.path(out_dir, "tables", "marker_rule_scores.tsv"))
write_tsv(marker_scores$cluster_labels, file.path(out_dir, "tables", "marker_rule_annotation.tsv"))

cluster_to_label <- setNames(marker_scores$cluster_labels$marker_rule_celltype, marker_scores$cluster_labels$cluster)
singlet$logic_marker_rule_celltype <- unname(cluster_to_label[as.character(singlet$seurat_clusters)])
singlet$logic_marker_rule_celltype[is.na(singlet$logic_marker_rule_celltype)] <- "Unknown"
singlet$logic_consensus_celltype <- singlet$logic_marker_rule_celltype

consensus <- data.frame(
  barcode = colnames(singlet),
  sample_id = singlet$sample_id,
  cluster = as.character(singlet$seurat_clusters),
  logic_marker_rule_celltype = singlet$logic_marker_rule_celltype,
  logic_consensus_celltype = singlet$logic_consensus_celltype,
  stringsAsFactors = FALSE
)
write_tsv(consensus, file.path(out_dir, "tables", "consensus_annotation.tsv"))

saveRDS(singlet, file.path(out_dir, "objects", "rds", "annotated_consensus.rds"))

p_cluster <- DimPlot(singlet, reduction = "umap", group.by = "seurat_clusters", label = TRUE, raster = FALSE) +
  ggtitle("Singlet clusters after doublet filtering") +
  theme_classic()
ggsave(file.path(out_dir, "figures", "umap_clusters_after_doublet.png"), p_cluster, width = 7, height = 5, dpi = 220)

p_annotation <- DimPlot(singlet, reduction = "umap", group.by = "logic_consensus_celltype", label = TRUE, raster = FALSE) +
  ggtitle("Marker-rule consensus annotation") +
  theme_classic()
ggsave(file.path(out_dir, "figures", "umap_marker_rule_annotation.png"), p_annotation, width = 8, height = 5.5, dpi = 220)

dot_features <- unique(unlist(marker_sets, use.names = FALSE))
dot_features <- intersect(dot_features, rownames(singlet))
if (length(dot_features) > 0) {
  p_dot <- DotPlot(singlet, features = dot_features, group.by = "logic_consensus_celltype") +
    RotatedAxis() +
    ggtitle("Consensus annotation marker evidence") +
    theme(axis.text.x = element_text(size = 8))
  ggsave(file.path(out_dir, "figures", "marker_dotplot_consensus.png"), p_dot, width = 12, height = 6, dpi = 220)
}

cell_counts <- data.frame(
  metric = c("cells_before_doublet_filter", "predicted_doublets", "cells_after_doublet_filter", "genes_after_filter", "clusters_after_filter"),
  value = c(ncol(obj), sum(predicted), ncol(singlet), nrow(singlet), length(levels(singlet$seurat_clusters)))
)
write_tsv(cell_counts, file.path(out_dir, "tables", "singlet_reprocess_summary.tsv"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "logs", "reprocess_singlets_sessionInfo.txt"))

append_status(out_dir, "post_doublet_reprocess", "PASS", paste(ncol(singlet), "singlets from", ncol(obj), "cells"))
append_status(out_dir, "marker_detection_singlet", "PASS", paste(nrow(markers), "marker rows"))
append_status(out_dir, "marker_rule_annotation", "PASS", paste(length(unique(singlet$logic_consensus_celltype)), "labels"))
cat("Reprocessed singlets:", nrow(singlet), "genes x", ncol(singlet), "cells\n")
