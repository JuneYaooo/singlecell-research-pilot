#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
})

MODULE <- "subset_extraction"

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

safe_name <- function(value) {
  cleaned <- gsub("[^A-Za-z0-9_.-]+", "_", value)
  cleaned <- gsub("^_+|_+$", "", cleaned)
  if (identical(cleaned, "")) "subset" else cleaned
}

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
out_dir <- args[["out"]]
if (is.null(out_dir)) {
  stop("Usage: run_logic_subset.R --out <logic_out> [--object <rds>] [--subset-column col] [--subset-values v1,v2] [--name name]")
}

input_rds <- value_or_default(args, "object", file.path(out_dir, "objects", "rds", "annotated_consensus.rds"))
subset_column <- value_or_default(args, "subset-column", "logic_consensus_celltype")
subset_values <- strsplit(value_or_default(args, "subset-values", "T_cells"), ",", fixed = TRUE)[[1]]
subset_values <- trimws(subset_values)
subset_values <- subset_values[nzchar(subset_values)]
subset_name <- safe_name(value_or_default(args, "name", paste(subset_values, collapse = "_")))

if (!file.exists(input_rds)) stop("Input RDS does not exist: ", input_rds)
if (length(subset_values) == 0) stop("--subset-values must contain at least one value")

dir.create(file.path(out_dir, "objects", "rds"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

append_status(out_dir, MODULE, "START", paste(subset_column, paste(subset_values, collapse = ","), sep = "="))

obj <- readRDS(input_rds)
if (!subset_column %in% colnames(obj@meta.data)) {
  stop("Subset column not found in Seurat metadata: ", subset_column)
}

meta <- obj@meta.data
keep <- rownames(meta)[as.character(meta[[subset_column]]) %in% subset_values]
if (length(keep) == 0) {
  stop("No cells matched ", subset_column, " in ", paste(subset_values, collapse = ","))
}

subset_obj <- subset(obj, cells = keep)
subset_obj@misc$logic_subset <- list(
  parent_object = input_rds,
  subset_column = subset_column,
  subset_values = subset_values,
  subset_name = subset_name,
  cells = length(keep)
)

object_path <- file.path(out_dir, "objects", "rds", paste0("subset_", subset_name, ".rds"))
saveRDS(subset_obj, object_path)

subset_meta <- subset_obj@meta.data
subset_meta$barcode <- rownames(subset_meta)
barcode_table <- data.frame(
  barcode = rownames(subset_meta),
  subset_name = subset_name,
  subset_column = subset_column,
  subset_value = as.character(subset_meta[[subset_column]]),
  stringsAsFactors = FALSE
)
write_tsv(barcode_table, file.path(out_dir, "tables", paste0("subset_", subset_name, "_barcodes.tsv")))

summary_columns <- intersect(c(subset_column, "sample_id", "condition", "seurat_clusters", "logic_consensus_celltype"), colnames(subset_meta))
summary <- subset_meta %>%
  mutate(across(all_of(summary_columns), as.character)) %>%
  group_by(across(all_of(summary_columns))) %>%
  summarise(cells = n(), .groups = "drop") %>%
  arrange(desc(cells))
write_tsv(summary, file.path(out_dir, "tables", paste0("subset_", subset_name, "_summary.tsv")))

params <- data.frame(
  parameter = c("parent_object", "subset_column", "subset_values", "subset_name", "cells", "output_object"),
  value = c(input_rds, subset_column, paste(subset_values, collapse = ","), subset_name, length(keep), object_path),
  stringsAsFactors = FALSE
)
write_tsv(params, file.path(out_dir, "tables", paste0("subset_", subset_name, "_params.tsv")))

detail <- paste0(
  "subset_name=", subset_name,
  " cells=", length(keep),
  " subset_column=", subset_column,
  " values=", paste(subset_values, collapse = ",")
)
append_status(out_dir, MODULE, "PASS", detail)
cat(detail, "\n")
