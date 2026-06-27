#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
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

args <- parse_args(commandArgs(trailingOnly = TRUE))
input_rds <- args[["input"]]
out_dir <- args[["out"]]

if (is.null(input_rds) || is.null(out_dir)) {
  stop("Usage: run_module_export_seurat.R --input <processed.rds> --out <output_root>")
}
if (!file.exists(input_rds)) {
  stop("Input RDS does not exist: ", input_rds)
}

dir.create(file.path(out_dir, "objects", "export"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "logs"), recursive = TRUE, showWarnings = FALSE)

append_status(out_dir, "export_seurat", "START", input_rds)

obj <- readRDS(input_rds)
DefaultAssay(obj) <- "RNA"
if (inherits(obj[["RNA"]], "Assay5")) {
  obj <- JoinLayers(obj, assay = "RNA")
}

counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
if (!inherits(counts, "sparseMatrix")) {
  counts <- as(counts, "dgCMatrix")
}

export_dir <- file.path(out_dir, "objects", "export")
plain_mtx <- file.path(export_dir, "counts.mtx")
gz_mtx <- paste0(plain_mtx, ".gz")
if (file.exists(plain_mtx)) unlink(plain_mtx)
if (file.exists(gz_mtx)) unlink(gz_mtx)
writeMM(counts, plain_mtx)
gzip_status <- system2("gzip", c("-f", plain_mtx), stdout = TRUE, stderr = TRUE)
if (!file.exists(gz_mtx)) {
  stop("Failed to gzip Matrix Market file: ", paste(gzip_status, collapse = " "))
}

write.table(
  data.frame(gene = rownames(counts)),
  file.path(export_dir, "genes.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
write.table(
  data.frame(barcode = colnames(counts)),
  file.path(export_dir, "barcodes.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

metadata <- obj@meta.data
metadata$barcode <- rownames(metadata)
metadata <- metadata[, c("barcode", setdiff(colnames(metadata), "barcode")), drop = FALSE]
write.table(
  metadata,
  file.path(export_dir, "metadata.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

clusters <- data.frame(
  barcode = colnames(obj),
  seurat_cluster = as.character(Idents(obj)),
  stringsAsFactors = FALSE
)
write.table(
  clusters,
  file.path(export_dir, "clusters.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

if ("umap" %in% Reductions(obj)) {
  umap <- as.data.frame(Embeddings(obj, "umap"))
  umap$barcode <- rownames(umap)
  umap <- umap[, c("barcode", setdiff(colnames(umap), "barcode")), drop = FALSE]
  write.table(
    umap,
    file.path(export_dir, "umap.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
} else {
  write.table(
    data.frame(barcode = colnames(obj)),
    file.path(export_dir, "umap.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

summary <- data.frame(
  metric = c("genes", "cells", "reductions", "clusters"),
  value = c(
    nrow(obj),
    ncol(obj),
    paste(Reductions(obj), collapse = ","),
    paste(levels(Idents(obj)), collapse = ",")
  )
)
write.table(
  summary,
  file.path(out_dir, "tables", "export_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "logs", "export_sessionInfo.txt"))
append_status(out_dir, "export_seurat", "PASS", paste(nrow(obj), "genes", ncol(obj), "cells"))
cat("Exported Seurat object:", nrow(obj), "genes x", ncol(obj), "cells\n")
