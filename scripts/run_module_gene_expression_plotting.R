#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(dplyr)
})

MODULE <- "gene_expression_plotting"

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
  stop("Usage: run_module_gene_expression_plotting.R --out <workflow_out> [--object <rds>] [--genes CD3D,MS4A1]")
}

input_rds <- value_or_default(args, "object", file.path(out_dir, "objects", "rds", "annotated_consensus.rds"))
genes_requested <- strsplit(value_or_default(args, "genes", "CD3D,MS4A1,EPCAM,LYZ,COL1A1,VWF"), ",", fixed = TRUE)[[1]]
genes_requested <- trimws(genes_requested)
genes_requested <- genes_requested[nzchar(genes_requested)]
group_by <- value_or_default(args, "group-by", "logic_consensus_celltype")
assay <- value_or_default(args, "assay", "RNA")
slot <- value_or_default(args, "slot", "data")
max_feature_genes <- as.integer(value_or_default(args, "max-feature-genes", "6"))

if (!file.exists(input_rds)) stop("Input RDS does not exist: ", input_rds)
if (length(genes_requested) == 0) stop("--genes must contain at least one gene")

dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)

append_status(out_dir, MODULE, "START", paste(genes_requested, collapse = ","))

obj <- readRDS(input_rds)
if (!assay %in% Assays(obj)) stop("Assay not found: ", assay)
DefaultAssay(obj) <- assay
if (inherits(obj[[assay]], "Assay5")) {
  obj <- JoinLayers(obj, assay = assay)
}
if (!group_by %in% colnames(obj@meta.data)) {
  stop("group-by column not found in metadata: ", group_by)
}
if (!"umap" %in% Reductions(obj)) {
  stop("Object does not contain a UMAP reduction")
}

feature_names <- rownames(obj[[assay]])
present <- intersect(genes_requested, feature_names)
missing <- setdiff(genes_requested, feature_names)
if (length(present) == 0) {
  stop("None of the requested genes are present in assay ", assay, ": ", paste(genes_requested, collapse = ","))
}

gene_status <- data.frame(
  gene = genes_requested,
  present = genes_requested %in% feature_names,
  assay = assay,
  slot = slot,
  stringsAsFactors = FALSE
)
write_tsv(gene_status, file.path(out_dir, "tables", "gene_expression_requested_genes.tsv"))

expr <- FetchData(obj, vars = present, layer = slot)
expr$barcode <- rownames(expr)
expr[[group_by]] <- as.character(obj@meta.data[rownames(expr), group_by])
summary_rows <- list()
for (gene in present) {
  stats <- expr %>%
    group_by(.data[[group_by]]) %>%
    summarise(
      gene = gene,
      cells = n(),
      mean_expression = mean(.data[[gene]], na.rm = TRUE),
      median_expression = median(.data[[gene]], na.rm = TRUE),
      pct_expressing = mean(.data[[gene]] > 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(group = all_of(group_by))
  summary_rows[[length(summary_rows) + 1]] <- stats
}
expression_summary <- bind_rows(summary_rows) %>%
  select(gene, group, cells, mean_expression, median_expression, pct_expressing)
write_tsv(expression_summary, file.path(out_dir, "tables", "gene_expression_summary.tsv"))

feature_genes <- head(present, max_feature_genes)
feature_plot <- FeaturePlot(
  obj,
  features = feature_genes,
  reduction = "umap",
  order = TRUE,
  ncol = min(3, length(feature_genes))
)
ggsave(
  file.path(out_dir, "figures", "gene_expression_featureplot.png"),
  feature_plot,
  width = max(5, min(12, 4 * min(3, length(feature_genes)))),
  height = max(4, 3.5 * ceiling(length(feature_genes) / min(3, length(feature_genes)))),
  dpi = 220
)

vln_plot <- VlnPlot(
  obj,
  features = feature_genes,
  group.by = group_by,
  pt.size = 0,
  ncol = min(3, length(feature_genes))
) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(
  file.path(out_dir, "figures", "gene_expression_vlnplot.png"),
  vln_plot,
  width = max(6, min(14, 4 * min(3, length(feature_genes)))),
  height = max(4, 3.8 * ceiling(length(feature_genes) / min(3, length(feature_genes)))),
  dpi = 220
)

dot_plot <- DotPlot(obj, features = present, group.by = group_by) +
  RotatedAxis() +
  labs(title = "Requested gene expression by annotated cell type") +
  theme_classic(base_size = 10)
ggsave(
  file.path(out_dir, "figures", "gene_expression_dotplot.png"),
  dot_plot,
  width = max(7, length(present) * 0.65),
  height = 4.8,
  dpi = 220
)

detail <- paste0(
  "requested_genes=", length(genes_requested),
  " present_genes=", length(present),
  " missing_genes=", length(missing),
  " group_by=", group_by,
  " assay=", assay,
  " slot=", slot
)
append_status(out_dir, MODULE, "PASS", detail)
cat(detail, "\n")
