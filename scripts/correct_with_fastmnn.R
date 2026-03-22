#!/usr/bin/env Rscript

# Suppress package startup messages for cleaner logs
suppressPackageStartupMessages({
  library(optparse)
  library(arrow)
  library(batchelor)
  library(SingleCellExperiment)
})

option_list = list(
  make_option(c("--input_data"), type = "character", help = "Path to input data file (parquet format)"),
  make_option(c("--batch_key"), type = "character", help = "Column name indicating batch information"),
  make_option(c("--output_path"), type = "character", help = "Path to save the corrected output file (parquet format)"),
  make_option(c("--k"), type = "integer", default = 20, help = "Number of nearest neighbors [default %default]"),
  make_option(c("--d"), type = "integer", default = 50, help = "Number of dimensions [default %default]"),
  make_option(c("--ndist"), type = "integer", default = 3, help = "Number of distances for variance estimation [default %default]"),
  make_option(c("--prop_k"), type = "double", default = NULL, help = "Proportion of cells for k [default NULL]"),
  make_option(c("--parameter_path"), type = "character", default = NULL, help = "Path to Optuna parameter CSV (overrides individual params)")
)

parser = OptionParser(option_list = option_list)
opt = parse_args(parser)

if (is.null(opt$input_data) || is.null(opt$batch_key) || is.null(opt$output_path)) {
  print_help(parser)
  stop("All three arguments (--input_data, --batch_key, --output_path) must be supplied.", call. = FALSE)
}

input_file = opt$input_data
batch_col = opt$batch_key
output_file = opt$output_path

# Load hyperparameters from Optuna CSV if provided
k_val = opt$k
d_val = opt$d
ndist_val = opt$ndist
prop_k_val = opt$prop_k

if (!is.null(opt$parameter_path)) {
  params_df = read.csv(opt$parameter_path)
  params_df = params_df[order(-params_df$total), ]
  best = params_df[1, ]
  if (best$state != "COMPLETE") stop("Optimization did not complete successfully")
  k_val = as.integer(best$params_k)
  d_val = as.integer(best$params_d)
  ndist_val = as.integer(best$params_ndist)
  prop_k_val = best$params_prop_k
  cat(sprintf("Using tuned params: k=%d, d=%d, ndist=%d, prop.k=%.3f\n", k_val, d_val, ndist_val, prop_k_val))
}

parquet_data = read_parquet(input_file)
col_names = names(parquet_data)
metadata_cols = col_names[grepl("^Metadata_", col_names)]
features_cols = col_names[!grepl("^Metadata_", col_names)]
features = parquet_data[, features_cols]
metadata = parquet_data[, metadata_cols]
batch_info = parquet_data[, batch_col]

mnn_args = list(t(features), batch=batch_info, k=k_val, d=d_val, ndist=ndist_val)
if (!is.null(prop_k_val)) {
  mnn_args$prop.k = prop_k_val
}
corrected = do.call(fastMNN, mnn_args)

corrected_df = cbind(metadata, reducedDim(corrected))
write_parquet(corrected_df, output_file)