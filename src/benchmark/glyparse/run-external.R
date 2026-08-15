#!/usr/bin/env Rscript

script_path <- sub(
  "^--file=",
  "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1L]]
)
source(file.path(dirname(normalizePath(script_path)), "common.R"))

benchmark_load_packages()
dir.create(benchmark_cache_dir, recursive = TRUE, showWarnings = FALSE)
force <- identical(Sys.getenv("GLYPARSE_BENCHMARK_FORCE", unset = "0"), "1")

run_python_tool <- function(tool, environment) {
  runtime_cache <- file.path(benchmark_cache_dir, paste0(tool, "-runtime"))
  dir.create(runtime_cache, recursive = TRUE, showWarnings = FALSE)
  arguments <- c(
    "run",
    "-n",
    environment,
    "python",
    file.path(benchmark_module_dir, "run-python-tool.py"),
    "--tool",
    tool,
    "--corpus-dir",
    benchmark_corpus_dir,
    "--output-dir",
    benchmark_cache_dir,
    "--workers",
    as.character(benchmark_cores),
    "--timeout",
    "5"
  )
  if (force) {
    arguments <- c(arguments, "--force")
  }
  status <- system2(
    "conda",
    arguments,
    env = c(
      paste0("MPLCONFIGDIR=", runtime_cache),
      paste0("XDG_CACHE_HOME=", runtime_cache)
    )
  )
  if (!identical(status, 0L)) {
    stop(tool, " conversion failed with exit status ", status)
  }
}

run_gfc <- function() {
  summary_path <- file.path(benchmark_cache_dir, "gfc-run-summary.csv")
  if (file.exists(summary_path) && !force) {
    message("Reusing completed GlycanFormatConverter raw conversion cache")
    return(invisible(NULL))
  }
  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    file.path(benchmark_module_dir, "run-gfc.R")
  )
  if (!identical(status, 0L)) {
    stop("GlycanFormatConverter conversion failed with exit status ", status)
  }
}

run_python_tool("glycowork", "glycowork")
run_python_tool("glypy", "glypy")
run_gfc()

tool_prefixes <- c(
  GlycanFormatConverter = "gfc",
  glycowork = "glycowork",
  glypy = "glypy"
)
formats <- setdiff(names(benchmark_definitions), "iupac_condensed")
raw_paths <- unlist(
  lapply(tool_prefixes, function(prefix) {
    file.path(benchmark_cache_dir, paste0(prefix, "-", formats, ".csv.gz"))
  }),
  use.names = FALSE
)
missing_paths <- raw_paths[!file.exists(raw_paths)]
if (length(missing_paths) > 0L) {
  stop(
    "Missing external result cache(s): ",
    paste(missing_paths, collapse = ", ")
  )
}

raw_outputs <- unlist(
  lapply(raw_paths, function(path) {
    data <- benchmark_read_csv(path)
    data$raw_iupac_condensed[data$conversion_status == "converted"]
  }),
  use.names = FALSE
)
raw_outputs <- raw_outputs[nzchar(raw_outputs)]
unique_outputs <- unique(raw_outputs)
message(
  "Normalizing ",
  length(unique_outputs),
  " unique external IUPAC outputs"
)

normalize_external_notation <- function(value) {
  value <- gsub("NeuAc", "Neu5Ac", value, fixed = TRUE)
  value <- gsub("NeuGc", "Neu5Gc", value, fixed = TRUE)
  value <- gsub("([[:alpha:]]+)2NAc", "\\1NAc", value, perl = TRUE)
  value
}

normalization_input <- vapply(
  unique_outputs,
  normalize_external_notation,
  character(1)
)
normalization_lookup <- data.frame(
  raw_iupac_condensed = unique_outputs,
  normalization_input = normalization_input,
  normalization_status = rep("", length(unique_outputs)),
  normalized_iupac_condensed = rep("", length(unique_outputs)),
  normalization_error = rep("", length(unique_outputs)),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
normalization_lookup_path <- file.path(
  benchmark_cache_dir,
  "external-normalization-lookup.csv.gz"
)
normalization_metadata_path <- file.path(
  benchmark_cache_dir,
  "external-normalization-metadata.json"
)
worktree_is_clean <- function(path) {
  output <- system2(
    "git",
    c("-C", path, "status", "--porcelain"),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(output, "status")
  (is.null(status) || status == 0L) && length(output) == 0L
}
normalization_revisions <- list(
  glyparse_revision = benchmark_revision(benchmark_glyparse_dir),
  glyrepr_revision = benchmark_revision(benchmark_glyrepr_dir),
  source_worktrees_clean = worktree_is_clean(benchmark_glyparse_dir) &&
    worktree_is_clean(benchmark_glyrepr_dir)
)
cache_hits <- rep(FALSE, length(unique_outputs))
if (
  file.exists(normalization_lookup_path) &&
    file.exists(normalization_metadata_path)
) {
  metadata <- jsonlite::read_json(normalization_metadata_path)
  revision_matches <- isTRUE(
    normalization_revisions$source_worktrees_clean
  ) && isTRUE(metadata$source_worktrees_clean) && identical(
    metadata$glyparse_revision,
    normalization_revisions$glyparse_revision
  ) && identical(
    metadata$glyrepr_revision,
    normalization_revisions$glyrepr_revision
  )
  if (revision_matches) {
    cached_lookup <- benchmark_read_csv(normalization_lookup_path)
    stopifnot(!anyDuplicated(cached_lookup$raw_iupac_condensed))
    cached_index <- match(
      normalization_lookup$raw_iupac_condensed,
      cached_lookup$raw_iupac_condensed
    )
    cache_hits <- !is.na(cached_index)
    cache_hits[cache_hits] <-
      cached_lookup$normalization_input[cached_index[cache_hits]] ==
        normalization_lookup$normalization_input[cache_hits]
    cached_rows <- which(cache_hits)
    cached_fields <- c(
      "normalization_status",
      "normalized_iupac_condensed",
      "normalization_error"
    )
    normalization_lookup[cached_rows, cached_fields] <-
      cached_lookup[cached_index[cached_rows], cached_fields]
  }
}
missing_rows <- which(!cache_hits)
message(
  "Reusing ",
  sum(cache_hits),
  " exact normalization result(s); computing ",
  length(missing_rows),
  " new result(s)"
)
if (length(missing_rows) > 0L) {
  normalized <- benchmark_normalize_iupac_fast(
    normalization_input[missing_rows]
  )
  normalization_lookup$normalization_status[missing_rows] <-
    normalized$normalization_status
  normalization_lookup$normalized_iupac_condensed[missing_rows] <-
    normalized$normalized_iupac_condensed
  normalization_lookup$normalization_error[missing_rows] <-
    normalized$normalization_error
}
stopifnot(all(nzchar(normalization_lookup$normalization_status)))
benchmark_write_csv(
  normalization_lookup,
  normalization_lookup_path
)
jsonlite::write_json(
  c(
    normalization_revisions,
    list(
      generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      rows = nrow(normalization_lookup),
      exact_cache_hits = sum(cache_hits),
      newly_normalized = length(missing_rows)
    )
  ),
  normalization_metadata_path,
  auto_unbox = TRUE,
  pretty = TRUE
)

summary_rows <- list()
summary_index <- 0L
for (tool in names(tool_prefixes)) {
  prefix <- tool_prefixes[[tool]]
  for (format in formats) {
    raw_path <- file.path(
      benchmark_cache_dir,
      paste0(prefix, "-", format, ".csv.gz")
    )
    result <- benchmark_read_csv(raw_path)
    match_index <- match(
      result$raw_iupac_condensed,
      normalization_lookup$raw_iupac_condensed
    )
    converted <- result$conversion_status == "converted"
    result$normalization_input <- ""
    result$normalization_status <- "not_attempted"
    result$normalized_iupac_condensed <- ""
    result$normalization_error <- ""
    available <- converted & !is.na(match_index)
    result$normalization_input[available] <-
      normalization_lookup$normalization_input[match_index[available]]
    result$normalization_status[available] <-
      normalization_lookup$normalization_status[match_index[available]]
    result$normalized_iupac_condensed[available] <-
      normalization_lookup$normalized_iupac_condensed[match_index[available]]
    result$normalization_error[available] <-
      normalization_lookup$normalization_error[match_index[available]]
    stopifnot(
      nrow(result) ==
        nrow(benchmark_read_corpus(benchmark_definitions[[format]])),
      all(!converted | available)
    )
    benchmark_write_csv(
      result,
      file.path(
        benchmark_cache_dir,
        paste0("external-", prefix, "-", format, ".csv.gz")
      )
    )
    summary_index <- summary_index + 1L
    summary_rows[[summary_index]] <- data.frame(
      tool = tool,
      tool_version = result$tool_version[[1L]],
      format = format,
      rows = nrow(result),
      converted = sum(converted),
      conversion_failed = sum(!converted),
      normalized = sum(result$normalization_status == "normalized"),
      normalization_failed = sum(result$normalization_status == "failed"),
      direct = sum(grepl("^direct", result$evidence_scope)),
      fallback = sum(
        result$evidence_scope == "accession_matched_wurcs_fallback"
      ),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
}

summary <- do.call(rbind, summary_rows)
benchmark_write_csv(
  summary,
  file.path(benchmark_cache_dir, "external-run-summary.csv")
)
print(summary)
