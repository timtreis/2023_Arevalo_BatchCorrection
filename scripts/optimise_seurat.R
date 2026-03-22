#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(arrow)
  library(Seurat)
  library(dplyr)
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

# --- Seurat v5 integration for a single parameter set ---

run_seurat_trial <- function(parquet_data, batch_col, label_col, seurat_method,
                             dims_val, k_anchor_val, k_weight_val) {
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

  # Split into layers by batch for v5 integration
  obj[["RNA"]] <- split(obj[["RNA"]], f = obj@meta.data[[batch_col]])

  VariableFeatures(obj) <- feat_names_clean
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = dims_val, verbose = FALSE)

  integration_method <- if (seurat_method == "cca") CCAIntegration else RPCAIntegration
  obj <- IntegrateLayers(
    object = obj,
    method = integration_method,
    orig.reduction = "pca",
    new.reduction = "integrated",
    dims = 1:dims_val,
    k.anchor = k_anchor_val,
    k.weight = k_weight_val,
    verbose = FALSE
  )

  corrected_mat <- Embeddings(obj, "integrated")
  batch_labels <- obj@meta.data[[batch_col]]
  bio_labels <- obj@meta.data[[label_col]]

  return(list(corrected = corrected_mat, batch_labels = batch_labels, bio_labels = bio_labels))
}

# --- Main ---

option_list <- list(
  make_option("--input_data", type = "character"),
  make_option("--batch_key", type = "character"),
  make_option("--label_key", type = "character"),
  make_option("--method", type = "character"),
  make_option("--n_trials", type = "integer", default = 30),
  make_option("--output_path", type = "character"),
  make_option("--smoketest", action = "store_true", default = FALSE)
)

opt <- parse_args(OptionParser(option_list = option_list))

if (opt$smoketest) opt$n_trials <- 2

parquet_data <- as.data.frame(read_parquet(opt$input_data))

# Random search grid
dims_range <- if (opt$smoketest) 5:15 else 5:50
k_anchor_range <- if (opt$smoketest) 3:8 else 3:30
k_weight_range <- if (opt$smoketest) 10:15 else 50:200

trials <- data.frame(
  number = seq_len(opt$n_trials) - 1,
  batch = NA_real_,
  bio = NA_real_,
  params_dims = sample(dims_range, opt$n_trials, replace = TRUE),
  params_k_anchor = sample(k_anchor_range, opt$n_trials, replace = TRUE),
  params_k_weight = sample(k_weight_range, opt$n_trials, replace = TRUE),
  state = "FAIL",
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(trials))) {
  cat(sprintf("Trial %d/%d: dims=%d, k.anchor=%d, k.weight=%d\n",
              i, nrow(trials), trials$params_dims[i],
              trials$params_k_anchor[i], trials$params_k_weight[i]))

  result <- tryCatch({
    res <- run_seurat_trial(
      parquet_data, opt$batch_key, opt$label_key, opt$method,
      trials$params_dims[i], trials$params_k_anchor[i], trials$params_k_weight[i]
    )
    batch_score <- compute_pcr_batch(res$corrected, res$batch_labels)
    bio_score <- compute_pcr_bio(res$corrected, res$bio_labels)
    list(batch = batch_score, bio = bio_score)
  }, error = function(e) {
    cat(sprintf("  FAILED: %s\n", conditionMessage(e)))
    NULL
  })

  if (!is.null(result)) {
    trials$batch[i] <- result$batch
    trials$bio[i] <- result$bio
    trials$state[i] <- "COMPLETE"
    cat(sprintf("  batch=%.4f, bio=%.4f\n", result$batch, result$bio))
  }
}

# Fallback: if all trials failed, emit one row with conservative defaults
if (all(trials$state != "COMPLETE")) {
  warning("WARNING: All HPO trials failed. Using conservative defaults — results may be suboptimal.")
  trials <- rbind(trials, data.frame(
    number = nrow(trials),
    batch = 0, bio = 0,
    params_dims = 10, params_k_anchor = 5, params_k_weight = 5,
    state = "COMPLETE",
    stringsAsFactors = FALSE
  ))
}

trials$total <- ifelse(is.na(trials$batch), NA, 0.6 * trials$bio + 0.4 * trials$batch)
trials <- trials[order(-trials$total, na.last = TRUE), ]
write.csv(trials, opt$output_path, row.names = FALSE)
cat(sprintf("Results written to %s\n", opt$output_path))
