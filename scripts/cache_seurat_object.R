#!/usr/bin/env Rscript
# Pre-builds a Seurat object from parquet and saves as .rds for fast trial loading.
# Performs: parquet load → matrix → Seurat object → split by batch → scale → PCA (max dims).

suppressPackageStartupMessages({
  library(optparse)
  library(arrow)
  library(Seurat)
  library(dplyr)
  library(future)
})

plan(sequential)
options(future.globals.maxSize = +Inf)

option_list <- list(
  make_option("--input_data", type = "character"),
  make_option("--batch_key", type = "character"),
  make_option("--max_dims", type = "integer", default = 50),
  make_option("--output_rds", type = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))

cat("Caching Seurat object from parquet...\n")

parquet_data <- as.data.frame(read_parquet(opt$input_data))
col_names <- names(parquet_data)
metadata_cols <- col_names[grepl("^Metadata_", col_names)]
features_cols <- col_names[!grepl("^Metadata_", col_names)]

feat_names_clean <- gsub("_", "-", features_cols)
expr_mat <- as.matrix(parquet_data[, features_cols])
colnames(expr_mat) <- feat_names_clean
rownames(expr_mat) <- paste0("cell_", seq_len(nrow(expr_mat)))
meta <- parquet_data[, metadata_cols, drop = FALSE]
rownames(meta) <- rownames(expr_mat)

# Free parquet_data early
rm(parquet_data)
gc()

obj <- CreateSeuratObject(counts = t(expr_mat), meta.data = meta)
obj <- SetAssayData(object = obj, layer = "data", new.data = t(expr_mat))

# Free expr_mat
rm(expr_mat)
gc()

obj[["RNA"]] <- split(obj[["RNA"]], f = obj@meta.data[[opt$batch_key]])
VariableFeatures(obj) <- feat_names_clean
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = opt$max_dims, verbose = FALSE)

saveRDS(obj, opt$output_rds)
cat(sprintf("Cached Seurat object to %s\n", opt$output_rds))
