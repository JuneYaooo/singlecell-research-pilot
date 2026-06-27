#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(harmony)
  library(hdf5r)
  library(R.utils)
})

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (length(hit) == 0 || hit == length(args)) return(default)
  args[[hit + 1]]
}

course_root_arg <- arg_value("--course-root", Sys.getenv("SEURAT_V5_COURSE_ROOT", unset = ""))
if (identical(course_root_arg, "")) {
  stop("--course-root or SEURAT_V5_COURSE_ROOT is required. Point it at the extracted Seurat V5 course directory.")
}
course_root <- normalizePath(course_root_arg, mustWork = TRUE)
out_root <- normalizePath(arg_value("--out", file.path(getwd(), "analysis", "seurat_v5_course_run")), mustWork = FALSE)
max_cells_per_sample <- as.integer(arg_value("--max-cells-per-sample", "2500"))
set.seed(20260626)

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "objects"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "logs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "tmp"), recursive = TRUE, showWarnings = FALSE)

status_path <- file.path(out_root, "tables", "run_status.tsv")
writeLines("step\tstatus\tdetail", status_path)
append_status <- function(step, status, detail = "") {
  detail <- gsub("[\r\n\t]+", " ", detail)
  cat(paste(step, status, detail, sep = "\t"), "\n", file = status_path, append = TRUE)
}

save_table <- function(x, name) {
  fwrite(as.data.table(x), file.path(out_root, "tables", name), sep = "\t")
}

save_plot <- function(plot, name, width = 8, height = 6) {
  ggsave(file.path(out_root, "figures", name), plot = plot, width = width, height = height)
}

add_qc <- function(obj, organism = "human") {
  mt_pattern <- if (organism == "mouse") "^mt-" else "^MT-"
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mt_pattern)
  obj[["percent.rb"]] <- PercentageFeatureSet(obj, pattern = "^RP[SL]")
  obj
}

object_summary <- function(obj, label) {
  data.frame(
    object = label,
    genes = nrow(obj),
    cells = ncol(obj),
    median_nFeature_RNA = median(obj$nFeature_RNA),
    median_nCount_RNA = median(obj$nCount_RNA),
    median_percent_mt = median(obj$percent.mt),
    stringsAsFactors = FALSE
  )
}

read_dense_count_table <- function(path, sep = "\t") {
  first_line <- readLines(gzfile(path), n = 1, warn = FALSE)
  first_fields <- strsplit(first_line, sep, fixed = TRUE)[[1]]
  header_has_gene_column <- length(first_fields) > 0 && identical(first_fields[[1]], "")
  if (header_has_gene_column) {
    dt <- fread(path, sep = sep, header = TRUE, data.table = TRUE, check.names = FALSE)
    gene_col <- names(dt)[1]
    genes <- dt[[gene_col]]
    dt[[gene_col]] <- NULL
  } else {
    dt <- fread(path, sep = sep, header = FALSE, skip = 1, data.table = TRUE, check.names = FALSE)
    genes <- dt[[1]]
    dt[[1]] <- NULL
    if (length(first_fields) != ncol(dt)) {
      stop("Header cell count (", length(first_fields), ") does not match matrix column count (", ncol(dt), ") for ", path)
    }
    setnames(dt, first_fields)
  }
  mat <- as.matrix(dt)
  mode(mat) <- "numeric"
  rownames(mat) <- make.unique(gsub("_", "-", as.character(genes)))
  colnames(mat) <- make.unique(as.character(colnames(mat)))
  Matrix(mat, sparse = TRUE)
}

normalize_10x_dir <- function(source_dir, target_dir) {
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  files <- list.files(source_dir, full.names = TRUE, recursive = TRUE)
  pick <- function(pattern) {
    hits <- files[grepl(pattern, basename(files), ignore.case = TRUE)]
    if (length(hits) == 0) stop("Missing 10X file matching ", pattern, " in ", source_dir)
    hits[[1]]
  }
  copy_as_10x <- function(src, plain_name) {
    dest_gz <- file.path(target_dir, paste0(plain_name, ".gz"))
    if (grepl("\\.gz$", src, ignore.case = TRUE)) {
      file.copy(src, dest_gz, overwrite = TRUE)
    } else {
      dest_plain <- file.path(target_dir, plain_name)
      file.copy(src, dest_plain, overwrite = TRUE)
      R.utils::gzip(dest_plain, destname = dest_gz, overwrite = TRUE, remove = TRUE)
    }
    dest_gz
  }
  old <- list.files(target_dir, full.names = TRUE)
  if (length(old) > 0) unlink(old)
  copy_as_10x(pick("(barcodes|count_matrix_barcodes).*tsv(\\.gz)?$"), "barcodes.tsv")
  copy_as_10x(pick("(features|genes|count_matrix_genes).*tsv(\\.gz)?$"), "features.tsv")
  copy_as_10x(pick("(matrix|sparse).*mtx(\\.gz)?$"), "matrix.mtx")
  target_dir
}

extract_tar_once <- function(tar_path, target_dir) {
  marker <- file.path(target_dir, ".extracted")
  if (!file.exists(marker)) {
    dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
    utils::untar(tar_path, exdir = target_dir)
    writeLines(as.character(Sys.time()), marker)
  }
  target_dir
}

run_step <- function(step, expr) {
  append_status(step, "START", "")
  tryCatch({
    value <- force(expr)
    append_status(step, "PASS", "")
    value
  }, error = function(e) {
    append_status(step, "FAIL", conditionMessage(e))
    stop(e)
  })
}

state <- new.env(parent = emptyenv())
state$all_summaries <- list()

single_obj <- run_step("03_single_matrix_import", {
  path <- file.path(course_root, "3.单矩阵读取", "GSE118389_counts_rsem.txt.gz")
  mat <- read_dense_count_table(path, sep = "\t")
  obj <- CreateSeuratObject(mat, project = "GSE118389_single_matrix", min.cells = 3, min.features = 0)
  obj <- add_qc(obj)
  state$all_summaries[["single_matrix"]] <- object_summary(obj, "single_matrix_GSE118389")
  saveRDS(obj, file.path(out_root, "objects", "single_matrix_import.rds"))
  obj
})

multi_matrix_objs <- run_step("04_multi_matrix_import", {
  path1 <- file.path(course_root, "4.多矩阵读取", "GSE118389_counts_rsem.txt.gz")
  path2 <- file.path(course_root, "4.多矩阵读取", "GSE158631_count.csv.gz")
  mat1 <- read_dense_count_table(path1, sep = "\t")
  mat2 <- read_dense_count_table(path2, sep = ",")
  obj1 <- add_qc(CreateSeuratObject(mat1, project = "GSE118389_multi_matrix", min.cells = 3, min.features = 0))
  obj2 <- add_qc(CreateSeuratObject(mat2, project = "GSE158631_multi_matrix", min.cells = 3, min.features = 0))
  state$all_summaries[["multi_matrix_1"]] <- object_summary(obj1, "multi_matrix_GSE118389")
  state$all_summaries[["multi_matrix_2"]] <- object_summary(obj2, "multi_matrix_GSE158631")
  saveRDS(obj1, file.path(out_root, "objects", "multi_matrix_GSE118389.rds"))
  saveRDS(obj2, file.path(out_root, "objects", "multi_matrix_GSE158631.rds"))
  list(obj1 = obj1, obj2 = obj2)
})

standard10x_objects <- run_step("05_standard_10x_tar_import", {
  tar_path <- file.path(course_root, "5.10X数据读取(标准10X)", "GSE184198_RAW.tar")
  extracted <- extract_tar_once(tar_path, file.path(out_root, "tmp", "GSE184198_RAW"))
  files <- list.files(extracted, pattern = "matrix.*mtx.gz$", full.names = TRUE)
  if (length(files) == 0) stop("No matrix files found after extracting ", tar_path)
  sample_prefixes <- unique(sub("(?:-|_)?matrix\\.mtx\\.gz$", "", basename(files)))
  objs <- lapply(sample_prefixes, function(prefix) {
    sample_files <- list.files(extracted, pattern = paste0("^", prefix), full.names = TRUE)
    src <- file.path(out_root, "tmp", paste0("standard10x_", make.names(prefix)))
    dir.create(src, recursive = TRUE, showWarnings = FALSE)
    file.copy(sample_files, src, overwrite = TRUE)
    norm <- normalize_10x_dir(src, file.path(out_root, "tmp", paste0("standard10x_norm_", make.names(prefix))))
    counts <- Read10X(norm)
    if (is.list(counts)) counts <- counts[[1]]
    obj <- CreateSeuratObject(counts, project = prefix, min.cells = 3, min.features = 200)
    obj$sample_id <- prefix
    obj$condition <- ifelse(grepl("GC", prefix), "GC", ifelse(grepl("NT", prefix), "NT", "unknown"))
    add_qc(obj)
  })
  names(objs) <- sample_prefixes
  state$all_summaries[["standard10x_tar"]] <- do.call(rbind, Map(object_summary, objs, names(objs)))
  saveRDS(objs, file.path(out_root, "objects", "standard10x_tar_import_list.rds"))
  objs
})

nonstandard10x_objects <- run_step("06_nonstandard_10x_tar_import", {
  tar_path <- file.path(course_root, "6.10X数据读取(非标准10X)", "GSE176078_RAW.tar")
  extracted <- extract_tar_once(tar_path, file.path(out_root, "tmp", "GSE176078_RAW"))
  inner <- list.files(extracted, pattern = "\\.tar\\.gz$", full.names = TRUE)
  if (length(inner) == 0) stop("No inner tar.gz files found after extracting ", tar_path)
  sample_dir <- extract_tar_once(inner[[1]], file.path(out_root, "tmp", tools::file_path_sans_ext(tools::file_path_sans_ext(basename(inner[[1]])))))
  norm <- normalize_10x_dir(sample_dir, file.path(out_root, "tmp", "GSE176078_first_norm"))
  counts <- Read10X(norm, gene.column = 1)
  if (is.list(counts)) counts <- counts[[1]]
  obj <- CreateSeuratObject(counts, project = "GSE176078_nonstandard_first", min.cells = 3, min.features = 200)
  obj$sample_id <- basename(inner[[1]])
  obj <- add_qc(obj)
  state$all_summaries[["nonstandard10x"]] <- object_summary(obj, "nonstandard10x_GSE176078_first")
  saveRDS(obj, file.path(out_root, "objects", "nonstandard10x_first_import.rds"))
  obj
})

multi10x_obj <- run_step("07_multiple_10x_import_merge", {
  sample_dirs <- c(
    file.path(course_root, "7.多个10X数据读取", "GSE184198", "GSM5580154"),
    file.path(course_root, "7.多个10X数据读取", "GSE184198", "GSM5580155")
  )
  objs <- lapply(sample_dirs, function(d) {
    counts <- Read10X(d)
    if (is.list(counts)) counts <- counts[[1]]
    sample_id <- basename(d)
    obj <- CreateSeuratObject(counts, project = sample_id, min.cells = 3, min.features = 200)
    obj$sample_id <- sample_id
    obj$condition <- ifelse(sample_id == "GSM5580154", "GC", "NT")
    obj <- add_qc(obj)
    cells <- colnames(obj)
    if (!is.na(max_cells_per_sample) && length(cells) > max_cells_per_sample) {
      obj <- subset(obj, cells = sample(cells, max_cells_per_sample))
    }
    obj
  })
  names(objs) <- basename(sample_dirs)
  merged <- merge(objs[[1]], y = objs[-1], add.cell.ids = names(objs), project = "GSE184198_multi10x")
  state$all_summaries[["multi10x_merged_raw"]] <- object_summary(merged, "multi10x_GSE184198_merged_downsampled")
  saveRDS(merged, file.path(out_root, "objects", "multi10x_merged_raw.rds"))
  merged
})

h5_obj <- run_step("08_h5_import", {
  tar_path <- file.path(course_root, "8.H5数据读取", "GSE203612_RAW.tar")
  extracted <- extract_tar_once(tar_path, file.path(out_root, "tmp", "GSE203612_RAW"))
  h5_files <- list.files(extracted, pattern = "\\.h5$", full.names = TRUE)
  if (length(h5_files) == 0) stop("No h5 files found after extracting ", tar_path)
  counts <- Read10X_h5(h5_files[[1]], use.names = TRUE)
  if (is.list(counts)) counts <- counts[[1]]
  obj <- CreateSeuratObject(counts, project = "GSE203612_h5_first", min.cells = 3, min.features = 200)
  obj$sample_id <- tools::file_path_sans_ext(basename(h5_files[[1]]))
  obj <- add_qc(obj)
  state$all_summaries[["h5_first"]] <- object_summary(obj, "h5_GSE203612_first")
  saveRDS(obj, file.path(out_root, "objects", "h5_first_import.rds"))
  obj
})

processed_obj <- run_step("10_11_13_qc_normalize_harmony_cluster", {
  obj <- multi10x_obj
  save_plot(VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4, pt.size = 0), "01_QC_before_violin.png", width = 12, height = 5)
  before <- ncol(obj)
  obj <- subset(obj, subset = nCount_RNA <= 25000 & nFeature_RNA >= 200 & nFeature_RNA <= 5000 & percent.mt <= 25 & percent.rb <= 40)
  after <- ncol(obj)
  save_table(data.frame(stage = c("before_qc", "after_qc"), cells = c(before, after)), "qc_cell_counts.tsv")
  save_plot(VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4, pt.size = 0), "02_QC_after_violin.png", width = 12, height = 5)

  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  save_plot(VariableFeaturePlot(obj), "03_variable_features.png", width = 7, height = 5)
  obj <- ScaleData(obj, features = rownames(obj), verbose = FALSE)
  obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 30, verbose = FALSE)
  save_plot(ElbowPlot(obj, ndims = 30), "04_elbow_plot.png", width = 6, height = 5)
  obj <- RunHarmony(obj, group.by.vars = "sample_id", reduction.use = "pca", dims.use = 1:20, verbose = FALSE)
  obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:20, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.2, verbose = FALSE)
  obj <- RunUMAP(obj, reduction = "harmony", dims = 1:20, verbose = FALSE)
  save_plot(DimPlot(obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE), "05_umap_clusters.png", width = 7, height = 6)
  save_plot(DimPlot(obj, reduction = "umap", group.by = "sample_id"), "06_umap_sample_id.png", width = 7, height = 6)
  save_plot(DimPlot(obj, reduction = "umap", group.by = "condition"), "07_umap_condition.png", width = 7, height = 6)
  state$all_summaries[["processed"]] <- object_summary(obj, "processed_multi10x_qc_harmony_cluster")
  saveRDS(obj, file.path(out_root, "objects", "processed_multi10x_qc_harmony_cluster.rds"))
  obj
})

markers <- run_step("17_marker_detection", {
  processed_for_markers <- JoinLayers(processed_obj)
  Idents(processed_for_markers) <- "seurat_clusters"
  markers <- FindAllMarkers(
    processed_for_markers,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25,
    test.use = "wilcox",
    max.cells.per.ident = 200,
    verbose = FALSE
  )
  fwrite(markers, file.path(out_root, "tables", "markers_FindAllMarkers.tsv"), sep = "\t")
  if (!"cluster" %in% colnames(markers)) {
    markers$cluster <- character()
  }
  top_markers <- markers
  if (nrow(markers) > 0) {
    top_markers <- markers |>
      dplyr::group_by(cluster) |>
      dplyr::slice_max(order_by = avg_log2FC, n = 5, with_ties = FALSE) |>
      dplyr::ungroup()
  }
  fwrite(top_markers, file.path(out_root, "tables", "markers_top5_by_cluster.tsv"), sep = "\t")
  if (nrow(top_markers) > 0) {
    genes <- unique(head(top_markers$gene, 20))
    save_plot(DotPlot(processed_for_markers, features = genes) + RotatedAxis(), "08_top_marker_dotplot.png", width = 10, height = 6)
  }
  markers
})

run_step("write_final_summary", {
  summary_df <- data.table::rbindlist(state$all_summaries, fill = TRUE)
  save_table(summary_df, "object_summaries.tsv")
  writeLines(capture.output(sessionInfo()), file.path(out_root, "logs", "sessionInfo.txt"))
  report <- c(
    "# Seurat V5 Course Core Run Report",
    "",
    paste("- Course root:", course_root),
    paste("- Output root:", out_root),
    paste("- Max cells per sample for main merged workflow:", max_cells_per_sample),
    paste("- Final processed object cells:", ncol(processed_obj)),
    paste("- Final processed object genes:", nrow(processed_obj)),
    paste("- Clusters:", paste(levels(processed_obj$seurat_clusters), collapse = ", ")),
    paste("- Marker rows:", nrow(markers)),
    "",
    "Validated course stages in this run:",
    "",
    "- 3.单矩阵读取",
    "- 4.多矩阵读取",
    "- 5.10X数据读取(标准10X)",
    "- 6.10X数据读取(非标准10X)",
    "- 7.多个10X数据读取",
    "- 8.H5数据读取",
    "- 10.数据质控",
    "- 11.标准化 / Harmony branch, excluding decontX because the core environment does not include celda/decontX",
    "- 13.分群聚类",
    "- 17.FindAllMarkers marker detection"
  )
  writeLines(report, file.path(out_root, "report.md"))
})

message("Completed Seurat V5 course core run: ", out_root)
