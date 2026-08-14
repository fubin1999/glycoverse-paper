#!/usr/bin/env Rscript

script_path <- sub(
  "^--file=",
  "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1L]]
)
source(file.path(dirname(normalizePath(script_path)), "common.R"))

gfc_protocol <- "glyparse-gfc-v1"
gfc_workers <- 4L
jar <- normalizePath(
  path.expand(Sys.getenv(
    "GLYCAN_FORMAT_CONVERTER_JAR",
    unset = "~/GlycanFormatConverter-cli/target/GlycanFormatConverter-cli.jar"
  )),
  mustWork = TRUE
)
reuse_path <- path.expand(Sys.getenv(
  "GLYPARSE_GFC_REUSE_CACHE",
  unset = paste0(
    "/Users/fubin/.codex/worktrees/cdc8/glyparse/docs/",
    "parser-validation/cache/gfc-reference-cache.rds"
  )
))
dir.create(benchmark_cache_dir, recursive = TRUE, showWarnings = FALSE)

gfc_encode <- function(value) {
  if (!nzchar(value)) {
    return("-")
  }
  jsonlite::base64_enc(charToRaw(enc2utf8(value)))
}

gfc_decode <- function(value) {
  if (identical(value, "-")) {
    return("")
  }
  rawToChar(jsonlite::base64_dec(value))
}

compile_runner <- function() {
  source_path <- file.path(benchmark_module_dir, "GfcBatchRunner.java")
  class_dir <- file.path(benchmark_cache_dir, "gfc-classes")
  class_path <- file.path(class_dir, "GfcBatchRunner.class")
  dir.create(class_dir, recursive = TRUE, showWarnings = FALSE)
  output <- system2(
    "javac",
    c(
      "-encoding",
      "UTF-8",
      "-cp",
      jar,
      "-d",
      class_dir,
      source_path
    ),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(output, "status")
  if ((!is.null(status) && status != 0L) || !file.exists(class_path)) {
    stop(
      "Unable to compile GfcBatchRunner.java:\n",
      paste(output, collapse = "\n")
    )
  }
  class_dir
}

write_requests <- function(requests, path) {
  lines <- paste(
    gfc_protocol,
    requests$request_id,
    vapply(requests$source_format, gfc_encode, character(1)),
    vapply(requests$source_sequence, gfc_encode, character(1)),
    sep = "\t"
  )
  writeLines(lines, path, useBytes = TRUE)
}

run_new_requests <- function(requests) {
  if (nrow(requests) == 0L) {
    return(data.frame())
  }
  class_dir <- compile_runner()
  class_path <- paste(class_dir, jar, sep = .Platform$path.sep)
  shard <- (requests$request_id - 1L) %% gfc_workers + 1L
  shard_indices <- split(seq_len(nrow(requests)), shard)
  run_shard <- function(index) {
    request_path <- file.path(
      benchmark_cache_dir,
      sprintf("gfc-requests-%02d.tsv", index)
    )
    output_path <- file.path(
      benchmark_cache_dir,
      sprintf("gfc-output-%02d.tsv", index)
    )
    if (file.exists(output_path)) {
      unlink(output_path)
    }
    write_requests(
      requests[shard_indices[[as.character(index)]], ],
      request_path
    )
    status <- system2(
      "java",
      c("-cp", class_path, "GfcBatchRunner", output_path),
      stdin = request_path,
      stdout = "",
      stderr = ""
    )
    if (!identical(status, 0L)) {
      stop("GfcBatchRunner shard ", index, " failed with status ", status)
    }
    output_path
  }
  active <- as.integer(names(shard_indices))
  output_paths <- parallel::mclapply(
    active,
    run_shard,
    mc.cores = min(length(active), gfc_workers),
    mc.preschedule = FALSE
  )
  if (any(vapply(output_paths, inherits, logical(1), "try-error"))) {
    stop("At least one GFC shard failed")
  }
  lines <- unlist(
    lapply(output_paths, readLines, warn = FALSE),
    use.names = FALSE
  )
  fields <- strsplit(lines, "\t", fixed = TRUE)
  stopifnot(all(lengths(fields) == 8L))
  result <- data.frame(
    request_id = as.integer(vapply(fields, `[[`, character(1), 2L)),
    conversion_status = vapply(fields, `[[`, character(1), 3L),
    failure_step = vapply(fields, `[[`, character(1), 4L),
    intermediate_wurcs = vapply(
      fields,
      function(x) gfc_decode(x[[5L]]),
      character(1)
    ),
    raw_iupac_condensed = vapply(
      fields,
      function(x) gfc_decode(x[[6L]]),
      character(1)
    ),
    conversion_warning = vapply(
      fields,
      function(x) gfc_decode(x[[7L]]),
      character(1)
    ),
    conversion_error = vapply(
      fields,
      function(x) gfc_decode(x[[8L]]),
      character(1)
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  result <- result[order(result$request_id), ]
  stopifnot(identical(result$request_id, requests$request_id))
  cbind(requests, result[setdiff(names(result), "request_id")])
}

reuse_gfc_cache <- function() {
  if (!file.exists(reuse_path)) {
    return(NULL)
  }
  cache <- readRDS(reuse_path)
  expected_hashes <- c(
    wurcs = benchmark_sha256(file.path(
      benchmark_corpus_dir,
      benchmark_definitions$wurcs$file
    )),
    glycoct = benchmark_sha256(file.path(
      benchmark_corpus_dir,
      benchmark_definitions$glycoct$file
    ))
  )
  actual_hashes <- c(
    wurcs = cache$metadata$corpora$wurcs$hash,
    glycoct = cache$metadata$corpora$glycoct$hash
  )
  if (
    !identical(expected_hashes, actual_hashes) ||
      !identical(benchmark_sha256(jar), cache$metadata$gfc$jar_hash)
  ) {
    stop(
      "The configured GFC reuse cache does not match the current corpora/jar"
    )
  }
  direct <- cache$results[
    cache$results$source_corpus %in% c("wurcs", "glycoct"),
  ]
  data.frame(
    format = direct$source_corpus,
    row_index = direct$source_row_index,
    glytoucan_ac = direct$glytoucan_ac,
    evidence_scope = "direct_source_converter_reused_verified_cache",
    conversion_status = direct$gfc_status,
    failure_step = direct$gfc_failure_step,
    intermediate_wurcs = direct$gfc_intermediate_wurcs,
    raw_iupac_condensed = direct$gfc_raw_iupac_condensed,
    conversion_warning = paste(
      direct$gfc_step1_stdout,
      direct$gfc_step1_stderr,
      direct$gfc_step2_stdout,
      direct$gfc_step2_stderr
    ),
    conversion_error = direct$gfc_runner_error,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

reused <- reuse_gfc_cache()
formats_to_run <- "iupac_extended"
if (is.null(reused)) {
  formats_to_run <- c("glycoct", "iupac_extended", "wurcs")
}
requests <- do.call(
  rbind,
  lapply(formats_to_run, function(format) {
    definition <- benchmark_definitions[[format]]
    data <- benchmark_read_corpus(definition)
    source_format <- switch(
      format,
      glycoct = "GlycoCT",
      iupac_extended = "IUPAC-Extended",
      wurcs = "WURCS"
    )
    data.frame(
      format = format,
      row_index = seq_len(nrow(data)),
      glytoucan_ac = data$glytoucan_ac,
      source_format = source_format,
      source_sequence = data[[definition$sequence_column]],
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
)
rownames(requests) <- NULL
requests$request_id <- seq_len(nrow(requests))
requests <- requests[c(
  "request_id",
  "format",
  "row_index",
  "glytoucan_ac",
  "source_format",
  "source_sequence"
)]

started <- proc.time()[["elapsed"]]
new_direct <- run_new_requests(requests)
elapsed <- proc.time()[["elapsed"]] - started
if (nrow(new_direct) > 0L) {
  new_direct$evidence_scope <- "direct_source_converter"
  new_direct <- new_direct[c(
    "format",
    "row_index",
    "glytoucan_ac",
    "evidence_scope",
    "conversion_status",
    "failure_step",
    "intermediate_wurcs",
    "raw_iupac_condensed",
    "conversion_warning",
    "conversion_error"
  )]
}
direct <- if (is.null(reused)) new_direct else rbind(reused, new_direct)
direct <- direct[order(direct$format, direct$row_index), ]

wurcs_direct <- direct[direct$format == "wurcs", ]
final_rows <- list()
for (format in setdiff(names(benchmark_definitions), "iupac_condensed")) {
  definition <- benchmark_definitions[[format]]
  corpus <- benchmark_read_corpus(definition)
  if (format %in% c("glycoct", "iupac_extended", "wurcs")) {
    result <- direct[direct$format == format, ]
  } else {
    matched <- match(corpus$glytoucan_ac, wurcs_direct$glytoucan_ac)
    result <- wurcs_direct[matched, ]
    result$format <- format
    result$row_index <- seq_len(nrow(corpus))
    result$glytoucan_ac <- corpus$glytoucan_ac
    result$evidence_scope <- "accession_matched_wurcs_fallback"
    missing <- is.na(matched)
    if (any(missing)) {
      result$conversion_status[missing] <- "unsupported"
      result$conversion_error[missing] <-
        "No accession-matched WURCS sequence is available."
    }
  }
  stopifnot(
    nrow(result) == nrow(corpus),
    identical(result$glytoucan_ac, corpus$glytoucan_ac)
  )
  result$tool <- "GlycanFormatConverter"
  result$tool_version <- "2.10.3"
  result <- result[c(
    "tool",
    "tool_version",
    "format",
    "row_index",
    "glytoucan_ac",
    "evidence_scope",
    "conversion_status",
    "raw_iupac_condensed",
    "conversion_warning",
    "conversion_error",
    "failure_step",
    "intermediate_wurcs"
  )]
  benchmark_write_csv(
    result,
    file.path(benchmark_cache_dir, paste0("gfc-", format, ".csv.gz"))
  )
  final_rows[[format]] <- data.frame(
    tool = "GlycanFormatConverter",
    tool_version = "2.10.3",
    format = format,
    rows = nrow(result),
    converted = sum(result$conversion_status == "converted"),
    failed = sum(result$conversion_status != "converted"),
    direct = sum(grepl("^direct", result$evidence_scope)),
    fallback = sum(result$evidence_scope == "accession_matched_wurcs_fallback"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

summary <- do.call(rbind, final_rows)
summary$new_conversion_elapsed_seconds <- elapsed
benchmark_write_csv(
  summary,
  file.path(benchmark_cache_dir, "gfc-run-summary.csv")
)
jsonlite::write_json(
  list(
    generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    jar = jar,
    jar_sha256 = benchmark_sha256(jar),
    version = benchmark_tool_version("java", c("-jar", jar, "-v")),
    reused_cache = if (is.null(reused)) "" else normalizePath(reuse_path),
    reused_cache_verified = !is.null(reused),
    new_conversion_elapsed_seconds = elapsed
  ),
  file.path(benchmark_cache_dir, "gfc-metadata.json"),
  auto_unbox = TRUE,
  pretty = TRUE
)
print(summary)
