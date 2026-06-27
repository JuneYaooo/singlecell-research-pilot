#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
})

MODULE <- "findmarkers_pairwise"

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

value_or_default <- function(values, key, default) {
  value <- values[[key]]
  if (is.null(value) || identical(value, "")) default else value
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

args <- parse_args(commandArgs(trailingOnly = TRUE))
out_dir <- args[["out"]]
if (is.null(out_dir)) {
  stop("Usage: run_module_findmarkers.R --out <workflow_out> [--object <rds>] [--group-by col] [--ident-1 A] [--ident-2 B]")
}

input_rds <- value_or_default(args, "object", file.path(out_dir, "objects", "rds", "annotated_consensus.rds"))
group_by <- value_or_default(args, "group-by", "logic_consensus_celltype")
ident_1 <- value_or_default(args, "ident-1", "T_cells")
ident_2 <- value_or_default(args, "ident-2", "B_cells")
assay <- value_or_default(args, "assay", "RNA")
min_pct <- as.numeric(value_or_default(args, "min-pct", "0.1"))
logfc_threshold <- as.numeric(value_or_default(args, "logfc-threshold", "0.25"))
test_use <- value_or_default(args, "test-use", "wilcox")
subset_column <- args[["subset-column"]]
subset_value <- args[["subset-value"]]

if (!file.exists(input_rds)) stop("Input RDS does not exist: ", input_rds)
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)

append_status(out_dir, MODULE, "START", paste(ident_1, "vs", ident_2, "by", group_by))

obj <- readRDS(input_rds)
if (!assay %in% Assays(obj)) stop("Assay not found: ", assay)
DefaultAssay(obj) <- assay
if (inherits(obj[[assay]], "Assay5")) {
  obj <- JoinLayers(obj, assay = assay)
}

if (!group_by %in% colnames(obj@meta.data)) {
  stop("group-by column not found in metadata: ", group_by)
}

if (!is.null(subset_column) && !is.null(subset_value)) {
  if (!subset_column %in% colnames(obj@meta.data)) {
    stop("subset-column not found in metadata: ", subset_column)
  }
  keep <- rownames(obj@meta.data)[as.character(obj@meta.data[[subset_column]]) == subset_value]
  if (length(keep) == 0) {
    stop("No cells matched subset ", subset_column, "=", subset_value)
  }
  obj <- subset(obj, cells = keep)
}

group_values <- as.character(obj@meta.data[[group_by]])
if (!ident_1 %in% group_values) stop("ident-1 not found in ", group_by, ": ", ident_1)
if (!ident_2 %in% group_values) stop("ident-2 not found in ", group_by, ": ", ident_2)

Idents(obj) <- group_values
cells_1 <- length(WhichCells(obj, idents = ident_1))
cells_2 <- length(WhichCells(obj, idents = ident_2))
if (cells_1 < 3 || cells_2 < 3) {
  stop("FindMarkers requires at least 3 cells per group; found ", ident_1, "=", cells_1, ", ", ident_2, "=", cells_2)
}

markers <- FindMarkers(
  obj,
  ident.1 = ident_1,
  ident.2 = ident_2,
  assay = assay,
  min.pct = min_pct,
  logfc.threshold = logfc_threshold,
  test.use = test_use
)
markers$gene <- rownames(markers)
markers <- markers[, c("gene", setdiff(colnames(markers), "gene")), drop = FALSE]
markers$comparison <- paste(ident_1, "vs", ident_2)
markers$group_by <- group_by
markers$cells_ident_1 <- cells_1
markers$cells_ident_2 <- cells_2
write_tsv(markers, file.path(out_dir, "tables", "findmarkers_pairwise.tsv"))

params <- data.frame(
  parameter = c(
    "object", "assay", "group_by", "ident_1", "ident_2", "cells_ident_1",
    "cells_ident_2", "min_pct", "logfc_threshold", "test_use",
    "subset_column", "subset_value"
  ),
  value = c(
    input_rds, assay, group_by, ident_1, ident_2, cells_1, cells_2,
    min_pct, logfc_threshold, test_use,
    ifelse(is.null(subset_column), "", subset_column),
    ifelse(is.null(subset_value), "", subset_value)
  ),
  stringsAsFactors = FALSE
)
write_tsv(params, file.path(out_dir, "tables", "findmarkers_pairwise_params.tsv"))

plot_df <- markers
plot_df$p_val_adj_plot <- pmax(as.numeric(plot_df$p_val_adj), .Machine$double.xmin)
plot_df$neg_log10_padj <- -log10(plot_df$p_val_adj_plot)
plot_df$direction <- ifelse(plot_df$avg_log2FC > 0, ident_1, ident_2)
plot_df$significant <- plot_df$p_val_adj < 0.05
plot_df$label <- ""
top_genes <- head(plot_df[order(plot_df$p_val_adj_plot, -abs(plot_df$avg_log2FC)), "gene"], 12)
plot_df$label[plot_df$gene %in% top_genes] <- plot_df$gene[plot_df$gene %in% top_genes]

volcano <- ggplot(plot_df, aes(x = avg_log2FC, y = neg_log10_padj)) +
  geom_point(aes(color = significant, shape = direction), alpha = 0.7, size = 1.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.3, color = "grey45") +
  geom_vline(xintercept = c(-logfc_threshold, logfc_threshold), linetype = "dotted", linewidth = 0.3, color = "grey45") +
  scale_color_manual(values = c("FALSE" = "#9e9e9e", "TRUE" = "#2c7fb8")) +
  labs(
    title = paste("FindMarkers:", ident_1, "vs", ident_2),
    x = paste0("avg_log2FC (+ favors ", ident_1, ")"),
    y = "-log10 adjusted p-value",
    color = "adj. p < 0.05",
    shape = "higher in"
  ) +
  theme_classic(base_size = 11)
ggsave(file.path(out_dir, "figures", "findmarkers_pairwise_volcano.png"), volcano, width = 7, height = 5, dpi = 220)

detail <- paste0(
  "comparison=", ident_1, "_vs_", ident_2,
  " group_by=", group_by,
  " rows=", nrow(markers),
  " cells_ident_1=", cells_1,
  " cells_ident_2=", cells_2
)
append_status(out_dir, MODULE, "PASS", detail)
cat(detail, "\n")
