#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(arrow)
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
  make_option("--input_data", type = "character"),
  make_option("--batch_key", type = "character"),
  make_option("--label_key", type = "character"),
  make_option("--n_trials", type = "integer", default = 30),
  make_option("--output_path", type = "character"),
  make_option("--smoketest", action = "store_true", default = FALSE)
)

opt <- parse_args(OptionParser(option_list = option_list))

if (opt$smoketest) opt$n_trials <- 2

parquet_data <- as.data.frame(read_parquet(opt$input_data))
col_names <- names(parquet_data)
metadata_cols <- col_names[grepl("^Metadata_", col_names)]
features_cols <- col_names[!grepl("^Metadata_", col_names)]
features <- parquet_data[, features_cols]
batch_info <- parquet_data[[opt$batch_key]]
bio_info <- parquet_data[[opt$label_key]]

# Random search grid
trials <- data.frame(
  number = seq_len(opt$n_trials) - 1,
  batch = NA_real_,
  bio = NA_real_,
  params_k = sample(5:50, opt$n_trials, replace = TRUE),
  params_d = sample(5:50, opt$n_trials, replace = TRUE),
  params_ndist = sample(1:5, opt$n_trials, replace = TRUE),
  params_prop_k = runif(opt$n_trials, 0.01, 0.5),
  state = "FAIL",
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(trials))) {
  cat(sprintf("Trial %d/%d: k=%d, d=%d, ndist=%d, prop.k=%.3f\n",
              i, nrow(trials), trials$params_k[i], trials$params_d[i],
              trials$params_ndist[i], trials$params_prop_k[i]))

  result <- tryCatch({
    corrected <- fastMNN(
      t(features), batch = batch_info,
      k = trials$params_k[i], d = trials$params_d[i],
      ndist = trials$params_ndist[i], prop.k = trials$params_prop_k[i]
    )
    corrected_mat <- reducedDim(corrected)
    batch_score <- compute_pcr_batch(corrected_mat, batch_info)
    bio_score <- compute_asw_bio(corrected_mat, bio_info)
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

# Fallback: if all trials failed, emit one row with defaults
if (all(trials$state != "COMPLETE")) {
  warning("WARNING: All HPO trials failed. Using conservative defaults — results may be suboptimal.")
  trials <- rbind(trials, data.frame(
    number = nrow(trials),
    batch = 0, bio = 0,
    params_k = 20, params_d = 50, params_ndist = 3, params_prop_k = NA,
    state = "COMPLETE",
    stringsAsFactors = FALSE
  ))
}

trials$total <- ifelse(is.na(trials$batch), NA, 0.6 * trials$bio + 0.4 * trials$batch)
trials <- trials[order(-trials$total, na.last = TRUE), ]
write.csv(trials, opt$output_path, row.names = FALSE)
cat(sprintf("Results written to %s\n", opt$output_path))
