#!/usr/bin/env Rscript
# Single-trial runner for fastMNN. Called by optimise_fastmnn.py via subprocess.
# Prints RESULT:<batch_score>,<bio_score> on success, exits 1 on failure.
# Supports --cached_rds for fast loading (pre-built data).

suppressPackageStartupMessages({
  library(optparse)
  library(batchelor)
  library(SingleCellExperiment)
})

# --- Lightweight evaluation metrics (pure R) ---

compute_pcr_batch <- function(corrected_mat, batch_labels) {
  n_pcs <- min(10, ncol(corrected_mat))
  pca <- prcomp(corrected_mat, center = TRUE, scale. = FALSE, rank. = n_pcs)
  total_var <- sum(pca$sdev^2)
  r2_sum <- 0
  for (i in seq_len(n_pcs)) {
    fit <- summary(lm(pca$x[, i] ~ as.factor(batch_labels)))
    r2_sum <- r2_sum + fit$r.squared * (pca$sdev[i]^2 / total_var)
  }
  return(1 - r2_sum)
}

compute_asw_bio <- function(corrected_mat, bio_labels) {
  n_pcs <- min(10, ncol(corrected_mat))
  pca <- prcomp(corrected_mat, center = TRUE, scale. = FALSE, rank. = n_pcs)
  total_var <- sum(pca$sdev^2)
  r2_sum <- 0
  for (i in seq_len(n_pcs)) {
    fit <- summary(lm(pca$x[, i] ~ as.factor(bio_labels)))
    r2_sum <- r2_sum + fit$r.squared * (pca$sdev[i]^2 / total_var)
  }
  return(r2_sum)
}

# --- Main ---

option_list <- list(
  make_option("--cached_rds", type = "character", default = NULL,
              help = "Path to cached .rds (skips parquet loading)"),
  make_option("--input_data", type = "character", default = NULL),
  make_option("--batch_key", type = "character"),
  make_option("--label_key", type = "character"),
  make_option("--k", type = "integer"),
  make_option("--d", type = "integer"),
  make_option("--ndist", type = "integer"),
  make_option("--prop_k", type = "double")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (!is.null(opt$cached_rds) && file.exists(opt$cached_rds)) {
  # Fast path: load pre-cached data
  cached <- readRDS(opt$cached_rds)
  features_t <- cached$features_t
  batch_info <- cached$batch_info
  bio_info <- cached$bio_info
} else {
  # Fallback: load from parquet
  suppressPackageStartupMessages(library(arrow))
  parquet_data <- as.data.frame(read_parquet(opt$input_data))
  col_names <- names(parquet_data)
  features_cols <- col_names[!grepl("^Metadata_", col_names)]
  features_t <- t(parquet_data[, features_cols])
  batch_info <- parquet_data[[opt$batch_key]]
  bio_info <- parquet_data[[opt$label_key]]
}

corrected <- fastMNN(
  features_t, batch = batch_info,
  k = opt$k, d = opt$d,
  ndist = opt$ndist, prop.k = opt$prop_k
)
corrected_mat <- reducedDim(corrected)
batch_score <- compute_pcr_batch(corrected_mat, batch_info)
bio_score <- compute_asw_bio(corrected_mat, bio_info)

cat(sprintf("RESULT:%.10f,%.10f\n", batch_score, bio_score))
