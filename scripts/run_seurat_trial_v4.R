#!/usr/bin/env Rscript
# Single-trial runner for Seurat v4 CCA/RPCA. Called by optimise_seurat_v4.py.
# Uses FindIntegrationAnchors + IntegrateData (Seurat v4 API).
# Prints RESULT:<batch_score>,<bio_score> on success, exits 1 on failure.

suppressPackageStartupMessages({
  library(optparse)
  library(arrow)
  library(Seurat)
})

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
  make_option("--input_data", type = "character"),
  make_option("--batch_key", type = "character"),
  make_option("--label_key", type = "character"),
  make_option("--method", type = "character"),
  make_option("--dims", type = "integer"),
  make_option("--k_anchor", type = "integer"),
  make_option("--k_weight", type = "integer")
)

opt <- parse_args(OptionParser(option_list = option_list))

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
  obj <- RunPCA(obj, features = colnames(raw), npcs = opt$dims, verbose = FALSE)
  seurat_lists[[i]] <- obj
}

result <- tryCatch({
  anchor_set <- FindIntegrationAnchors(
    object.list = seurat_lists,
    reduction = opt$method,
    anchor.features = gsub("_", "-", features_cols),
    dims = 1:opt$dims,
    k.anchor = opt$k_anchor,
    scale = FALSE,
    verbose = FALSE
  )
  IntegrateData(
    anchorset = anchor_set,
    dims = 1:opt$dims,
    k.weight = opt$k_weight,
    verbose = FALSE
  )
}, error = function(e) {
  cat(sprintf("Seurat v4 integration failed: %s\n", conditionMessage(e)), file = stderr())
  NULL
})

if (is.null(result)) {
  quit(status = 1)
}

int_data <- GetAssayData(object = result, assay = "integrated", slot = "data")
corrected_mat <- int_data %>% as.matrix() %>% t()

batch_labels <- result@meta.data[[opt$batch_key]]
bio_labels <- result@meta.data[[opt$label_key]]

batch_score <- compute_pcr_batch(corrected_mat, batch_labels)
bio_score <- compute_pcr_bio(corrected_mat, bio_labels)

cat(sprintf("RESULT:%.10f,%.10f\n", batch_score, bio_score))
