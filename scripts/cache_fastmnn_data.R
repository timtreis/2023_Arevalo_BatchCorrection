#!/usr/bin/env Rscript
# Pre-loads parquet data and saves as .rds for fast trial loading.
# Caches: transposed feature matrix, batch vector, bio labels vector.

suppressPackageStartupMessages({
  library(optparse)
  library(arrow)
})

option_list <- list(
  make_option("--input_data", type = "character"),
  make_option("--batch_key", type = "character"),
  make_option("--label_key", type = "character"),
  make_option("--output_rds", type = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))

cat("Caching fastMNN data from parquet...\n")

parquet_data <- as.data.frame(read_parquet(opt$input_data))
col_names <- names(parquet_data)
features_cols <- col_names[!grepl("^Metadata_", col_names)]

cached <- list(
  features_t = t(parquet_data[, features_cols]),
  batch_info = parquet_data[[opt$batch_key]],
  bio_info = parquet_data[[opt$label_key]]
)

saveRDS(cached, opt$output_rds)
cat(sprintf("Cached fastMNN data to %s\n", opt$output_rds))
