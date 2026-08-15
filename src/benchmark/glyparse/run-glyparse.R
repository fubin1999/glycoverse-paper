#!/usr/bin/env Rscript

script_path <- sub(
  "^--file=",
  "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1L]]
)
source(file.path(dirname(normalizePath(script_path)), "common.R"))

benchmark_load_packages()
dir.create(benchmark_cache_dir, recursive = TRUE, showWarnings = FALSE)

reference_definition <- benchmark_definitions$iupac_condensed
reference <- benchmark_read_corpus(reference_definition)
reference_normalized <- benchmark_normalize_iupac(
  reference[[reference_definition$sequence_column]]
)
reference_table <- data.frame(
  glytoucan_ac = reference$glytoucan_ac,
  reference_iupac_condensed = reference[[reference_definition$sequence_column]],
  reference_normalization_status = reference_normalized$normalization_status,
  reference_normalized_iupac_condensed = reference_normalized$normalized_iupac_condensed,
  reference_normalization_error = reference_normalized$normalization_error,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
benchmark_write_csv(
  reference_table,
  file.path(benchmark_cache_dir, "reference-normalized.csv.gz")
)

benchmark_rows <- list()
run_rows <- list()

for (format in names(benchmark_definitions)) {
  definition <- benchmark_definitions[[format]]
  data <- benchmark_read_corpus(definition)
  parsed_run <- benchmark_parse_corpus(format, definition, data)
  result <- cbind(
    data.frame(
      format = format,
      format_label = definition$label,
      parser = definition$parser,
      row_index = seq_len(nrow(data)),
      glytoucan_ac = data$glytoucan_ac,
      source_sequence = data[[definition$sequence_column]],
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    parsed_run$result
  )

  match_index <- match(result$glytoucan_ac, reference_table$glytoucan_ac)
  has_reference <- !is.na(match_index)
  result$reference_validation_applicable <-
    definition$reference_validation & has_reference
  result$reference_iupac_condensed <- ""
  result$reference_normalization_status <- "not_applicable"
  result$reference_normalized_iupac_condensed <- ""
  result$reference_normalization_error <- ""
  if (definition$reference_validation && any(has_reference)) {
    result$reference_iupac_condensed[has_reference] <-
      reference_table$reference_iupac_condensed[match_index[has_reference]]
    result$reference_normalization_status[has_reference] <-
      reference_table$reference_normalization_status[match_index[has_reference]]
    result$reference_normalized_iupac_condensed[has_reference] <-
      reference_table$reference_normalized_iupac_condensed[match_index[
        has_reference
      ]]
    result$reference_normalization_error[has_reference] <-
      reference_table$reference_normalization_error[match_index[has_reference]]
  }
  comparable <- result$parse_status == "parsed" &
    result$reference_validation_applicable &
    result$reference_normalization_status == "normalized"
  result$reference_comparison_status <- if (!definition$reference_validation) {
    rep("not_requested", nrow(result))
  } else {
    equivalent <- rep(FALSE, nrow(result))
    equivalent[comparable] <- benchmark_semantically_equal(
      result$parsed_iupac_condensed[comparable],
      result$reference_normalized_iupac_condensed[comparable]
    )
    ifelse(
      !has_reference,
      "no_accession_reference",
      ifelse(
        result$reference_normalization_status != "normalized",
        "reference_normalization_failed",
        ifelse(
          result$parse_status != "parsed",
          "parser_failed",
          ifelse(
            equivalent,
            "equivalent",
            "semantic_difference"
          )
        )
      )
    )
  }
  stopifnot(
    nrow(result) == nrow(data),
    identical(result$glytoucan_ac, data$glytoucan_ac),
    all(
      !comparable |
        result$reference_comparison_status %in%
          c("equivalent", "semantic_difference")
    )
  )

  benchmark_write_csv(
    result,
    file.path(benchmark_cache_dir, paste0("glyparse-", format, ".csv.gz"))
  )
  benchmark_rows[[format]] <- benchmark_measure_parser(
    format,
    definition,
    data,
    parsed_run$result
  )
  run_rows[[format]] <- data.frame(
    format = format,
    rows = nrow(result),
    parsed = sum(result$parse_status == "parsed"),
    failed = sum(result$parse_status == "failed"),
    validation_elapsed_seconds = parsed_run$elapsed,
    validation_rows_per_second = nrow(result) / parsed_run$elapsed,
    reference_available = sum(result$reference_validation_applicable),
    reference_equivalent = sum(
      result$reference_comparison_status == "equivalent"
    ),
    reference_semantic_difference = sum(
      result$reference_comparison_status == "semantic_difference"
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

benchmark_results <- do.call(rbind, benchmark_rows)
rownames(benchmark_results) <- NULL
run_summary <- do.call(rbind, run_rows)
rownames(run_summary) <- NULL
benchmark_write_csv(
  benchmark_results,
  file.path(benchmark_cache_dir, "glyparse-benchmark-repeats.csv")
)
benchmark_write_csv(
  run_summary,
  file.path(benchmark_cache_dir, "glyparse-run-summary.csv")
)

metadata <- list(
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  r_version = R.version.string,
  glyparse_version = as.character(utils::packageVersion("glyparse")),
  glyparse_revision = benchmark_revision(benchmark_glyparse_dir),
  glyrepr_version = as.character(utils::packageVersion("glyrepr")),
  glyrepr_revision = benchmark_revision(benchmark_glyrepr_dir),
  corpus = lapply(benchmark_definitions, function(definition) {
    path <- file.path(benchmark_corpus_dir, definition$file)
    list(
      file = definition$file,
      sha256 = benchmark_sha256(path),
      rows = nrow(benchmark_read_corpus(definition))
    )
  })
)
jsonlite::write_json(
  metadata,
  file.path(benchmark_cache_dir, "glyparse-metadata.json"),
  auto_unbox = TRUE,
  pretty = TRUE
)

print(run_summary)
