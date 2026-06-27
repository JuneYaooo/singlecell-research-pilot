#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})

MODULE <- "export_consensus_embeddings"

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
  stop("Usage: run_logic_export_consensus_embeddings.R --out <logic_out> [--object <rds>]")
}

input_rds <- value_or_default(args, "object", file.path(out_dir, "objects", "rds", "annotated_consensus.rds"))
if (!file.exists(input_rds)) stop("Input RDS does not exist: ", input_rds)

dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
append_status(out_dir, MODULE, "START", input_rds)

obj <- readRDS(input_rds)

metadata <- obj@meta.data
metadata$barcode <- rownames(metadata)
metadata <- metadata[, c("barcode", setdiff(colnames(metadata), "barcode")), drop = FALSE]
metadata$active_ident <- as.character(Idents(obj))
write_tsv(metadata, file.path(out_dir, "tables", "annotated_consensus_metadata.tsv"))

if (!"umap" %in% Reductions(obj)) {
  stop("annotated consensus object does not contain a umap reduction")
}
umap <- as.data.frame(Embeddings(obj, "umap"))
if (ncol(umap) < 2) stop("UMAP reduction has fewer than 2 dimensions")
umap <- umap[, 1:2, drop = FALSE]
colnames(umap) <- c("umap_1", "umap_2")
umap$barcode <- rownames(umap)
umap <- umap[, c("barcode", "umap_1", "umap_2"), drop = FALSE]
write_tsv(umap, file.path(out_dir, "tables", "annotated_consensus_umap.tsv"))

reductions <- names(obj@reductions)
reduction_rows <- data.frame(
  reduction = reductions,
  dimensions = vapply(reductions, function(name) ncol(Embeddings(obj, name)), integer(1)),
  cells = vapply(reductions, function(name) nrow(Embeddings(obj, name)), integer(1)),
  stringsAsFactors = FALSE
)
write_tsv(reduction_rows, file.path(out_dir, "tables", "annotated_consensus_reductions.tsv"))

detail <- paste0(
  "cells=", ncol(obj),
  " genes=", nrow(obj),
  " metadata_cols=", ncol(metadata),
  " reductions=", paste(reductions, collapse = ",")
)
append_status(out_dir, MODULE, "PASS", detail)
cat(detail, "\n")
