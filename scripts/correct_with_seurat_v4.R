#!/usr/bin/env Rscript
# Seurat v4 correction using FindIntegrationAnchors + IntegrateData (Arevalo et al. API).

suppressPackageStartupMessages({
  library(optparse)
  library(arrow)
  library(Seurat)
  library(dplyr)
})

option_list <- list(
  make_option(c("--input_data"), type = "character", help = "Path to input data file"),
  make_option(c("--batch_key"), type = "character", help = "Batch column name"),
  make_option(c("--method"), type = "character", help = "Seurat integration method (cca or rpca)"),
  make_option(c("--output_path"), type = "character", help = "Path to output file"),
  make_option(c("--dims"), type = "integer", default = 30, help = "Number of dimensions [default %default]"),
  make_option(c("--k_anchor"), type = "integer", default = 5, help = "Neighbors for anchor finding [default %default]"),
  make_option(c("--k_weight"), type = "integer", default = 100, help = "Neighbors for weighting [default %default]"),
  make_option(c("--parameter_path"), type = "character", default = NULL, help = "Path to HPO parameter CSV")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)
if (is.null(opt$input_data) || is.null(opt$batch_key) || is.null(opt$method) || is.null(opt$output_path)) {
  print_help(parser)
  stop("Required arguments missing.", call. = FALSE)
}

dims_val <- opt$dims
k_anchor_val <- opt$k_anchor
k_weight_val <- opt$k_weight

if (!is.null(opt$parameter_path)) {
  params_df <- read.csv(opt$parameter_path)
  complete_df <- params_df[params_df$state == "COMPLETE", ]
  if (nrow(complete_df) == 0) stop("No COMPLETE trials in HPO -- fix upstream before correction")
  complete_df <- complete_df[order(-complete_df$total), ]
  best <- complete_df[1, ]
  dims_val <- as.integer(best$params_dims)
  k_anchor_val <- as.integer(best$params_k_anchor)
  k_weight_val <- as.integer(best$params_k_weight)
  cat(sprintf("Using tuned params: dims=%d, k.anchor=%d, k.weight=%d\n", dims_val, k_anchor_val, k_weight_val))
}

options(future.globals.maxSize = +Inf)

parquet_data <- as.data.frame(read_parquet(opt$input_data))
col_names <- names(parquet_data)
metadata_cols <- col_names[grepl("^Metadata_", col_names)]
features_cols <- col_names[!grepl("^Metadata_", col_names)]

batch_info <- parquet_data[, opt$batch_key]
batches <- split(parquet_data[, features_cols], batch_info)
meta_batches <- split(parquet_data[, metadata_cols], batch_info)

seurat_lists <- list()
for (i in seq_along(batches)) {
  raw <- as.data.frame(batches[[i]])
  meta <- as.data.frame(meta_batches[[i]])

  names(raw) <- gsub("_", "-", names(raw))
  obj <- CreateSeuratObject(counts = t(raw), meta.data = meta)
  obj <- SetAssayData(object = obj, slot = "data", new.data = t(as.matrix(raw)))
  obj <- SetAssayData(object = obj, slot = "scale.data", new.data = t(as.matrix(raw)))
  obj <- RunPCA(obj, features = colnames(raw), npcs = dims_val, verbose = FALSE)

  seurat_lists[[i]] <- obj
}

anchor_set <- FindIntegrationAnchors(
  object.list = seurat_lists,
  reduction = opt$method,
  anchor.features = gsub("_", "-", features_cols),
  dims = 1:dims_val,
  k.anchor = k_anchor_val,
  scale = FALSE,
  verbose = FALSE
)

integrated_obj <- IntegrateData(
  anchorset = anchor_set,
  dims = 1:dims_val,
  k.weight = k_weight_val,
  verbose = FALSE
)

int_data <- GetAssayData(object = integrated_obj, assay = "integrated", slot = "data")
corrected <- int_data %>% as.matrix() %>% t() %>% as.data.frame()
colnames(corrected) <- gsub("-", "_", colnames(corrected))
meta_corrected <- integrated_obj@meta.data[, metadata_cols, drop = FALSE]

result_df <- cbind(meta_corrected, corrected)
write_parquet(result_df, opt$output_path)
cat(sprintf("Seurat v4 %s integration complete. Output: %s\n", opt$method, opt$output_path))
