#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(Matrix)
  library(limma)
  library(e1071)
  library(preprocessCore)
  library(reshape2)
  library(ggplot2)
  library(doParallel)
})

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (length(hit) == 0 || hit == length(args)) return(default)
  args[[hit + 1]]
}

course_dir_arg <- arg_value("--course-dir", "")
if (identical(course_dir_arg, "")) {
  course_root_arg <- Sys.getenv("SINGLECELL_COURSE_ROOT", unset = "")
  if (identical(course_root_arg, "")) {
    stop("--course-dir or SINGLECELL_COURSE_ROOT is required. Point it at the CIBERSORT course folder or extracted Seurat course directory.")
  }
  course_dir_arg <- file.path(course_root_arg, "24-1反卷积(CIBERSORT)")
}
course_dir <- normalizePath(course_dir_arg, mustWork = TRUE)
seurat_rds <- normalizePath(arg_value("--seurat-rds", file.path(
  getwd(), "analysis", "course_run", "objects", "processed_multi10x_qc_harmony_cluster.rds"
)), mustWork = TRUE)
out_dir <- normalizePath(arg_value("--out", file.path(getwd(), "analysis", "course_run", "cibersort")), mustWork = FALSE)
max_samples <- as.integer(arg_value("--max-samples", "30"))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
status_path <- file.path(out_dir, "cibersort_status.tsv")
writeLines("step\tstatus\tdetail", status_path)
status <- function(step, state, detail = "") {
  detail <- gsub("[\r\n\t]+", " ", detail)
  cat(paste(step, state, detail, sep = "\t"), "\n", file = status_path, append = TRUE)
}

run_step <- function(step, expr) {
  status(step, "START")
  tryCatch({
    value <- force(expr)
    status(step, "PASS")
    value
  }, error = function(e) {
    status(step, "FAIL", conditionMessage(e))
    stop(e)
  })
}

obj <- run_step("load_processed_seurat", {
  readRDS(seurat_rds)
})

signature <- run_step("build_cluster_signature_from_seurat", {
  obj <- JoinLayers(obj)
  counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
  clusters <- as.character(obj$seurat_clusters)
  sig <- sapply(sort(unique(clusters)), function(cl) {
    cells <- colnames(obj)[clusters == cl]
    Matrix::rowMeans(counts[, cells, drop = FALSE])
  })
  colnames(sig) <- paste0("cluster_", colnames(sig))
  sig <- as.data.frame(sig)
  sig$ID <- rownames(sig)
  setDT(sig)
  setcolorder(sig, "ID")
  lmgenes <- fread(file.path(course_dir, "LMgenes.txt"))
  sig <- sig[ID %in% lmgenes$ID]
  if (nrow(sig) < 50) stop("Too few intersecting LMgenes for CIBERSORT signature: ", nrow(sig))
  fwrite(sig, file.path(out_dir, "scimmune_from_clusters.txt"), sep = "\t", quote = FALSE)
  sig
})

mixture <- run_step("prepare_tcga_mixture_subset", {
  tcga <- fread(file.path(course_dir, "TCGAexp.txt"), data.table = FALSE, check.names = FALSE)
  rownames(tcga) <- tcga[[1]]
  tcga <- tcga[, -1, drop = FALSE]
  tumor_cols <- grep("-01A", colnames(tcga), value = TRUE, fixed = TRUE)
  if (length(tumor_cols) == 0) stop("No TCGA tumor sample columns ending in -01A found")
  tumor_cols <- head(tumor_cols, max_samples)
  sub <- tcga[, tumor_cols, drop = FALSE]
  colnames(sub) <- sub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1-\\2-\\3", colnames(sub))
  out <- rbind(ID = colnames(sub), sub)
  fwrite(as.data.table(out, keep.rownames = TRUE), file.path(out_dir, "uniq.symbol.subset.txt"), sep = "\t", quote = FALSE, col.names = FALSE)
  data.frame(genes = nrow(sub), samples = ncol(sub))
})

results <- run_step("run_cibersort", {
  old <- getwd()
  setwd(out_dir)
  on.exit(setwd(old), add = TRUE)
  Sys.setenv(OBJC_DISABLE_INITIALIZE_FORK_SAFETY = "YES")
  detectCores <- function(...) 1L
  source(file.path(course_dir, "xxdimmune.R"), encoding = "utf-8")
  CIBERSORT("scimmune_from_clusters.txt", "uniq.symbol.subset.txt", perm = 2, QN = FALSE)
  res <- fread("CIBERSORT-Results.txt")
  fwrite(res, "CIBERSORT-Results.tsv", sep = "\t")
  res
})

run_step("plot_cibersort_overview", {
  fraction_cols <- setdiff(colnames(results), c("Mixture", "P-value", "Correlation", "RMSE"))
  frac <- as.data.frame(results[, ..fraction_cols])
  rownames(frac) <- results$Mixture
  frac <- as.matrix(frac)
  png(file.path(out_dir, "cibersort_fraction_barplot.png"), width = 1600, height = 900, res = 140)
  par(las = 2, mar = c(10, 5, 4, 12))
  barplot(t(frac), col = rainbow(ncol(frac)), ylab = "Relative fraction", cex.names = 0.65)
  legend("topright", legend = colnames(frac), fill = rainbow(ncol(frac)), cex = 0.75, xpd = TRUE, inset = c(-0.18, 0))
  dev.off()
})

run_step("write_cibersort_report", {
  report <- c(
    "# CIBERSORT Course Case Run",
    "",
    paste("- Course dir:", course_dir),
    paste("- Seurat object:", seurat_rds),
    paste("- Signature LMgene rows:", nrow(signature)),
    paste("- TCGA genes:", mixture$genes),
    paste("- TCGA subset samples:", mixture$samples),
    paste("- Result rows:", nrow(results)),
    paste("- Result columns:", ncol(results))
  )
  writeLines(report, file.path(out_dir, "report.md"))
  writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
})

message("Completed CIBERSORT course case run: ", out_dir)
