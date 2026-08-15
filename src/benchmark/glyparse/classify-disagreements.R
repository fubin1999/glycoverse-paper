#!/usr/bin/env Rscript

script_path <- sub(
  "^--file=",
  "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1L]]
)
source(file.path(dirname(normalizePath(script_path)), "common.R"))

path <- file.path(benchmark_output_dir, "semantic_disagreements.csv.gz")
disagreements <- benchmark_read_csv(path)
classified <- benchmark_classify_semantic_differences(disagreements)

stopifnot(
  nrow(classified) == nrow(disagreements),
  all(nzchar(classified$primary_mismatch_class)),
  all(classified$mismatch_class_count >= 1L)
)

benchmark_write_csv(classified, path)
message("Classified ", nrow(classified), " semantic disagreement rows.")
