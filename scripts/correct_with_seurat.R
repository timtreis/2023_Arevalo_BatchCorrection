#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(arrow)
  library(Seurat)
  library(dplyr)
  library(future)
})

option_list <- list(
  make_option(c("--input_data"), type = "character", help = "Path to input data file"),
  make_option(c("--batch_key"), type = "character", help = "Batch column name"),
  make_option(c("--method"), type = "character", help = "Seurat integration method"),
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
  params_df <- params_df[order(-params_df$total), ]
  best <- params_df[1, ]
  if (best$state != "COMPLETE") stop("Optimization did not complete successfully")
  dims_val <- as.integer(best$params_dims)
  k_anchor_val <- as.integer(best$params_k_anchor)
  k_weight_val <- as.integer(best$params_k_weight)
  cat(sprintf("Using tuned params: dims=%d, k.anchor=%d, k.weight=%d\n", dims_val, k_anchor_val, k_weight_val))
}

# Force sequential plan to avoid serializing large globals to parallel workers
plan(sequential)
options(future.globals.maxSize = +Inf)

parquet_data <- as.data.frame(read_parquet(opt$input_data))
col_names <- names(parquet_data)
metadata_cols <- col_names[grepl("^Metadata_", col_names)]
features_cols <- col_names[!grepl("^Metadata_", col_names)]

# Clean feature names for R compatibility
feat_names_clean <- gsub("_", "-", features_cols)

# Build a single Seurat object with batch as a metadata column, then split into layers
expr_mat <- as.matrix(parquet_data[, features_cols])
colnames(expr_mat) <- feat_names_clean
rownames(expr_mat) <- paste0("cell_", seq_len(nrow(expr_mat)))
meta <- parquet_data[, metadata_cols, drop = FALSE]
rownames(meta) <- rownames(expr_mat)

obj <- CreateSeuratObject(counts = t(expr_mat), meta.data = meta)
obj <- SetAssayData(object = obj, layer = "data", new.data = t(expr_mat))

# Seurat v5: split into layers by batch for integration
obj[["RNA"]] <- split(obj[["RNA"]], f = obj@meta.data[[opt$batch_key]])

# Standard workflow: find variable features, scale, PCA
VariableFeatures(obj) <- feat_names_clean
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = dims_val, verbose = FALSE)

# Seurat v5 integration via IntegrateLayers
obj <- IntegrateLayers(
  object = obj,
  method = if (opt$method == "cca") CCAIntegration else RPCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated",
  dims = 1:dims_val,
  k.anchor = k_anchor_val,
  k.weight = k_weight_val,
  verbose = FALSE
)

# Extract integrated embedding
integrated_embed <- Embeddings(obj, "integrated")
corrected_df <- as.data.frame(integrated_embed)
colnames(corrected_df) <- paste0("seurat_", seq_len(ncol(corrected_df)))

# Rejoin layers and combine with metadata
meta_out <- obj@meta.data[, metadata_cols, drop = FALSE]
result_df <- cbind(meta_out, corrected_df)

write_parquet(result_df, opt$output_path)
cat(sprintf("Seurat %s integration complete. Output: %s\n", opt$method, opt$output_path))
