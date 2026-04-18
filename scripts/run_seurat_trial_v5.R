#!/usr/bin/env Rscript
# Single-trial runner for Seurat CCA/RPCA. Called by optimise_seurat.py via subprocess.
# Prints RESULT:<batch_score>,<bio_score> on success, exits 1 on failure.
# Supports --cached_rds for fast loading (pre-built Seurat object with PCA).

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(future)
})

# Force sequential plan to avoid serializing large globals to parallel workers
plan(sequential)
options(future.globals.maxSize = +Inf)

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

compute_pcr_bio <- function(corrected_mat, bio_labels) {
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
              help = "Path to cached Seurat .rds (skips parquet loading)"),
  make_option("--input_data", type = "character", default = NULL),
  make_option("--batch_key", type = "character"),
  make_option("--label_key", type = "character"),
  make_option("--method", type = "character"),
  make_option("--dims", type = "integer"),
  make_option("--k_anchor", type = "integer"),
  make_option("--k_weight", type = "integer")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (!is.null(opt$cached_rds) && file.exists(opt$cached_rds)) {
  # Fast path: load pre-built Seurat object
  obj <- readRDS(opt$cached_rds)
} else {
  # Fallback: build from parquet (backwards compatible)
  suppressPackageStartupMessages(library(arrow))
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

  obj <- CreateSeuratObject(counts = t(expr_mat), meta.data = meta)
  obj <- SetAssayData(object = obj, layer = "data", new.data = t(expr_mat))
  obj[["RNA"]] <- split(obj[["RNA"]], f = obj@meta.data[[opt$batch_key]])
  VariableFeatures(obj) <- feat_names_clean
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = opt$dims, verbose = FALSE)
}

integration_method <- if (opt$method == "cca") CCAIntegration else RPCAIntegration

# Retry integration with progressively lower k.weight if it fails or produces NaN.
# This handles cross-microscope scenarios where anchor counts vary across source pairs.
k_weight_try <- opt$k_weight
integrated_obj <- NULL

while (k_weight_try >= 3) {
  result <- tryCatch(
    IntegrateLayers(
      object = obj,
      method = integration_method,
      orig.reduction = "pca",
      new.reduction = "integrated",
      dims = 1:opt$dims,
      k.anchor = opt$k_anchor,
      k.weight = k_weight_try,
      verbose = FALSE
    ),
    error = function(e) {
      cat(sprintf("IntegrateLayers failed with k.weight=%d: %s\n",
                  k_weight_try, conditionMessage(e)), file = stderr())
      NULL
    }
  )

  if (!is.null(result)) {
    mat <- Embeddings(result, "integrated")
    if (all(is.finite(mat))) {
      integrated_obj <- result
      break
    } else {
      n_bad <- sum(!is.finite(mat))
      cat(sprintf("k.weight=%d produced %d non-finite values, retrying lower\n",
                  k_weight_try, n_bad), file = stderr())
    }
  }

  k_weight_try <- k_weight_try %/% 2
}

if (is.null(integrated_obj)) {
  cat("All k.weight attempts failed\n", file = stderr())
  quit(status = 1)
}

if (k_weight_try != opt$k_weight) {
  cat(sprintf("EFFECTIVE_K_WEIGHT:%d\n", k_weight_try))
}

corrected_mat <- Embeddings(integrated_obj, "integrated")
batch_labels <- integrated_obj@meta.data[[opt$batch_key]]
bio_labels <- integrated_obj@meta.data[[opt$label_key]]

batch_score <- compute_pcr_batch(corrected_mat, batch_labels)
bio_score <- compute_pcr_bio(corrected_mat, bio_labels)

cat(sprintf("RESULT:%.10f,%.10f\n", batch_score, bio_score))
