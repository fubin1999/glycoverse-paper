#!/usr/bin/env Rscript

script_path <- sub(
  "^--file=",
  "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1L]]
)
source(file.path(dirname(normalizePath(script_path)), "common.R"))

benchmark_load_packages()
dir.create(benchmark_output_dir, recursive = TRUE, showWarnings = FALSE)

external_tools <- c(
  GlycanFormatConverter = "gfc",
  glycowork = "glycowork",
  glypy = "glypy"
)
formats <- names(benchmark_definitions)
validated_formats <- setdiff(formats, "iupac_condensed")

collapse_counts <- function(values) {
  values <- values[nzchar(values)]
  if (length(values) == 0L) {
    return("none")
  }
  counts <- sort(table(values))
  paste0(names(counts), " x", as.integer(counts), collapse = "; ")
}

canonical_features <- function(value) {
  floating <- regmatches(
    value,
    gregexpr("\\{[^{}]*\\}", value, perl = TRUE)
  )[[1L]]
  residues <- regmatches(
    value,
    gregexpr("[[:alnum:]?+-]+(?=\\()", value, perl = TRUE)
  )[[1L]]
  linkages <- regmatches(
    value,
    gregexpr("(?<=\\()[ab?][0-9?/]+-[0-9?/]*", value, perl = TRUE)
  )[[1L]]
  floating_domains <- if (identical(floating, "")) {
    character()
  } else {
    ifelse(
      grepl("\\|[^|{}]+\\}$", floating),
      sub("^.*\\|([^|{}]+)\\}$", "\\1", floating),
      "all-main-nodes"
    )
  }
  list(
    residues = if (identical(residues, "")) character() else residues,
    linkages = if (identical(linkages, "")) character() else linkages,
    floating = if (identical(floating, "")) character() else floating,
    floating_domains = floating_domains,
    alditol = grepl("-ol\\(", value),
    reducing_linkage = if (length(linkages) == 0L) "" else tail(linkages, 1L),
    branch_shape = gsub(
      "[[:alnum:]?+-]+\\([ab?][0-9?/]+-[0-9?/]*\\)?",
      "R",
      gsub("\\{[^{}]*\\}", "F", value, perl = TRUE),
      perl = TRUE
    )
  )
}

canonical_difference_reason <- function(left, right) {
  if (!nzchar(left) || !nzchar(right)) {
    return(
      "One side has no normalized glycan, so no semantic comparison is possible."
    )
  }
  if (identical(left, right)) {
    return("The normalized canonical glycans are identical.")
  }
  a <- canonical_features(left)
  b <- canonical_features(right)
  reasons <- character()
  if (!identical(sort(a$residues), sort(b$residues))) {
    reasons <- c(
      reasons,
      paste0(
        "Residue or residue-substituent composition differs (glyparse: ",
        collapse_counts(a$residues),
        "; comparator: ",
        collapse_counts(b$residues),
        ")."
      )
    )
  }
  if (!identical(sort(a$linkages), sort(b$linkages))) {
    reasons <- c(
      reasons,
      paste0(
        "Linkage or anomer multiset differs (glyparse: ",
        collapse_counts(a$linkages),
        "; comparator: ",
        collapse_counts(b$linkages),
        ")."
      )
    )
  }
  if (!identical(a$alditol, b$alditol)) {
    reasons <- c(reasons, "Reducing-end alditol state differs.")
  }
  if (!identical(a$reducing_linkage, b$reducing_linkage)) {
    reasons <- c(
      reasons,
      paste0(
        "Reducing-end anomer/linkage differs (glyparse: ",
        a$reducing_linkage,
        "; comparator: ",
        b$reducing_linkage,
        ")."
      )
    )
  }
  if (!identical(sort(a$floating), sort(b$floating))) {
    reasons <- c(
      reasons,
      paste0(
        "Floating component/substituent or candidate-parent domain differs ",
        "(glyparse: ",
        length(a$floating),
        " part(s), domains ",
        collapse_counts(a$floating_domains),
        "; comparator: ",
        length(b$floating),
        " part(s), domains ",
        collapse_counts(b$floating_domains),
        ")."
      )
    )
  }
  if (length(reasons) == 0L || !identical(a$branch_shape, b$branch_shape)) {
    reasons <- c(
      reasons,
      "Canonical topology or residue-to-residue attachment differs after node-order normalization."
    )
  }
  paste(unique(reasons), collapse = " ")
}

diagnose_normalized_source <- function(format, source_sequence) {
  graph_parser <- switch(
    format,
    glycoct = function(value) {
      glyparse:::do_parse_glycoct(gsub(" ", "\n", value, fixed = TRUE))
    },
    wurcs = glyparse:::do_parse_wurcs,
    NULL
  )
  if (!is.null(graph_parser)) {
    graph <- tryCatch(graph_parser(source_sequence), error = identity)
    if (inherits(graph, "condition")) {
      return(c("source_parse", benchmark_condition_text(graph)))
    }
    validation <- tryCatch(
      glyrepr:::validate_glycan_graph(graph),
      error = identity
    )
    if (inherits(validation, "condition")) {
      return(c("glyrepr_validation", benchmark_condition_text(validation)))
    }
    return(c(
      "wrapper_validation",
      "The internal source parse and graph validation succeeded during replay; the exported parser failure is therefore wrapper- or context-dependent."
    ))
  }
  normalizer <- switch(
    format,
    gwb = glyparse:::convert_one_gwb_to_condensed,
    iupac_compact = glyparse:::convert_iupac_compact_to_condensed,
    iupac_condensed = identity,
    iupac_extended = function(value) {
      glyparse:::convert_ext_to_con(
        glyparse:::normalize_iupac_extended(value)
      )
    },
    NULL
  )
  if (is.null(normalizer)) {
    return(c(
      "",
      "No current internal normalizer is available for diagnostic replay."
    ))
  }
  normalized <- tryCatch(normalizer(source_sequence), error = identity)
  if (inherits(normalized, "condition")) {
    return(c("normalization", benchmark_condition_text(normalized)))
  }
  result <- tryCatch(
    glyrepr::as_glycan_structure(normalized, on_failure = "error"),
    error = identity
  )
  if (inherits(result, "condition")) {
    return(c(
      "glyrepr_construction",
      paste0(
        "Normalized to ",
        encodeString(normalized, quote = '"'),
        "; glyrepr construction failed: ",
        benchmark_condition_text(result)
      )
    ))
  }
  c(
    "wrapper_validation",
    "The internal conversion and glyrepr construction succeeded during replay; the exported parser failure is therefore wrapper- or context-dependent."
  )
}

classify_failure <- function(reason, historical_class = "") {
  if (nzchar(historical_class)) {
    return(historical_class)
  }
  if (
    grepl(
      paste(
        "must be either all generic or all concrete",
        "out tree",
        "cannot be represented",
        "glyrepr",
        "Invalid linkage",
        "no conflict-free",
        "repeating-unit",
        sep = "|"
      ),
      reason,
      ignore.case = TRUE
    )
  ) {
    return("glyrepr_restriction")
  }
  "glyparse_limitation"
}

historical_dir <- Sys.getenv(
  "GLYPARSE_PRIOR_DIAGNOSTIC_DIR",
  unset = paste0(
    "/Users/fubin/.codex/worktrees/cdc8/glyparse/",
    "docs/parser-validation"
  )
)

historical_failure_lookup <- function(format, current) {
  path <- file.path(historical_dir, paste0("parse_", format, ".csv"))
  if (!file.exists(path)) {
    return(NULL)
  }
  historical <- benchmark_read_csv(path)
  key <- paste(current$glytoucan_ac, current$source_sequence, sep = "\r")
  historical_key <- paste(
    historical$glytoucan_ac,
    historical$source_sequence,
    sep = "\r"
  )
  index <- match(key, historical_key)
  data.frame(
    failure_stage = ifelse(is.na(index), "", historical$failure_stage[index]),
    failure_class = ifelse(is.na(index), "", historical$failure_class[index]),
    failure_reason = ifelse(is.na(index), "", historical$failure_reason[index]),
    diagnostic_source = ifelse(
      is.na(index) | !nzchar(historical$failure_reason[index]),
      "",
      "exact-row prior diagnostic replay"
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

diagnose_failures <- function(format, failures) {
  detailed <- historical_failure_lookup(format, failures)
  if (is.null(detailed)) {
    detailed <- data.frame(
      failure_stage = rep("", nrow(failures)),
      failure_class = rep("", nrow(failures)),
      failure_reason = rep("", nrow(failures)),
      diagnostic_source = rep("", nrow(failures)),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
  needs_replay <- !nzchar(detailed$failure_reason)
  if (any(needs_replay)) {
    replay <- parallel::mclapply(
      failures$source_sequence[needs_replay],
      function(value) diagnose_normalized_source(format, value),
      mc.cores = benchmark_cores,
      mc.preschedule = TRUE
    )
    replay <- do.call(rbind, replay)
    detailed$failure_stage[needs_replay] <- replay[, 1L]
    detailed$failure_reason[needs_replay] <- replay[, 2L]
    detailed$diagnostic_source[
      needs_replay
    ] <- "current internal diagnostic replay"
  }
  detailed$failure_class <- vapply(
    seq_len(nrow(detailed)),
    function(index) {
      classify_failure(
        detailed$failure_reason[[index]],
        detailed$failure_class[[index]]
      )
    },
    character(1)
  )
  detailed
}

row_results <- list()
failure_results <- list()
difference_results <- list()
format_summaries <- list()
tool_summaries <- list()
difference_index <- 0L

for (format in formats) {
  message("Adjudicating ", format)
  parser <- benchmark_read_csv(file.path(
    benchmark_cache_dir,
    paste0("glyparse-", format, ".csv.gz")
  ))
  n <- nrow(parser)
  row <- parser[, c(
    "format",
    "format_label",
    "parser",
    "row_index",
    "glytoucan_ac",
    "source_sequence",
    "parse_status",
    "parsed_iupac_condensed",
    "parse_error",
    "reference_validation_applicable",
    "reference_normalized_iupac_condensed",
    "reference_comparison_status"
  )]

  reference_differs <- which(
    parser$reference_comparison_status == "semantic_difference"
  )
  if (length(reference_differs) > 0L) {
    difference_index <- difference_index + 1L
    difference_results[[difference_index]] <- data.frame(
      format = format,
      row_index = parser$row_index[reference_differs],
      glytoucan_ac = parser$glytoucan_ac[reference_differs],
      source_sequence = parser$source_sequence[reference_differs],
      comparator = "IUPAC-condensed accession reference",
      evidence_scope = "accession_matched_reference",
      glyparse_iupac_condensed = parser$parsed_iupac_condensed[
        reference_differs
      ],
      comparator_raw_iupac_condensed = parser$reference_iupac_condensed[
        reference_differs
      ],
      comparator_normalized_iupac_condensed = parser$reference_normalized_iupac_condensed[
        reference_differs
      ],
      difference_reason = vapply(
        reference_differs,
        function(index) {
          canonical_difference_reason(
            parser$parsed_iupac_condensed[[index]],
            parser$reference_normalized_iupac_condensed[[index]]
          )
        },
        character(1)
      ),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  external_canonical <- list()
  external_scope <- list()
  for (tool in names(external_tools)) {
    prefix <- external_tools[[tool]]
    if (identical(format, "iupac_condensed")) {
      next
    }
    external <- benchmark_read_csv(file.path(
      benchmark_cache_dir,
      paste0("external-", prefix, "-", format, ".csv.gz")
    ))
    stopifnot(
      nrow(external) == n,
      identical(external$row_index, parser$row_index),
      identical(external$glytoucan_ac, parser$glytoucan_ac)
    )
    canonical <- ifelse(
      external$normalization_status == "normalized",
      external$normalized_iupac_condensed,
      ""
    )
    external_canonical[[tool]] <- canonical
    external_scope[[tool]] <- external$evidence_scope
    key <- gsub("[^[:alnum:]]", "_", tolower(tool))
    row[[paste0(key, "_conversion_status")]] <- external$conversion_status
    row[[paste0(key, "_normalization_status")]] <-
      external$normalization_status
    row[[paste0(key, "_normalized_iupac_condensed")]] <- canonical
    row[[paste0(key, "_evidence_scope")]] <- external$evidence_scope
    comparable <- parser$parse_status == "parsed" & nzchar(canonical)
    comparison <- rep("not_comparable", n)
    comparison[comparable] <- ifelse(
      parser$parsed_iupac_condensed[comparable] == canonical[comparable],
      "equivalent",
      "semantic_difference"
    )
    row[[paste0(key, "_comparison_status")]] <- comparison

    differs <- which(comparison == "semantic_difference")
    if (length(differs) > 0L) {
      difference_index <- difference_index + 1L
      difference_results[[difference_index]] <- data.frame(
        format = format,
        row_index = parser$row_index[differs],
        glytoucan_ac = parser$glytoucan_ac[differs],
        source_sequence = parser$source_sequence[differs],
        comparator = tool,
        evidence_scope = external$evidence_scope[differs],
        glyparse_iupac_condensed = parser$parsed_iupac_condensed[differs],
        comparator_raw_iupac_condensed = external$raw_iupac_condensed[differs],
        comparator_normalized_iupac_condensed = canonical[differs],
        difference_reason = vapply(
          differs,
          function(index) {
            canonical_difference_reason(
              parser$parsed_iupac_condensed[[index]],
              canonical[[index]]
            )
          },
          character(1)
        ),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
    tool_summaries[[paste(format, tool)]] <- data.frame(
      format = format,
      tool = tool,
      rows = n,
      converted = sum(external$conversion_status %in% c("converted")),
      conversion_failed = sum(!external$conversion_status %in% c("converted")),
      normalized = sum(external$normalization_status == "normalized"),
      normalization_failed = sum(external$normalization_status == "failed"),
      comparable_with_glyparse = sum(comparable),
      equivalent_to_glyparse = sum(comparison == "equivalent"),
      different_from_glyparse = sum(comparison == "semantic_difference"),
      direct_rows = sum(grepl("^direct", external$evidence_scope)),
      fallback_rows = sum(
        external$evidence_scope == "accession_matched_wurcs_fallback"
      ),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  row$adjudication_status <- if (identical(format, "iupac_condensed")) {
    ifelse(row$parse_status == "parsed", "parse_coverage_only", "parser_failed")
  } else {
    rep("parsed_unvalidated", n)
  }
  row$adjudication_basis <- if (identical(format, "iupac_condensed")) {
    ifelse(
      row$parse_status == "parsed",
      "Correctness validation was not requested for IUPAC-condensed source sequences.",
      "The parser did not produce a glycan."
    )
  } else {
    rep(
      "No normalized accession reference or external result was available.",
      n
    )
  }

  if (!identical(format, "iupac_condensed")) {
    external_matrix <- do.call(cbind, external_canonical)
    scope_matrix <- do.call(cbind, external_scope)
    parser_failed <- row$parse_status == "failed"
    row$adjudication_status[parser_failed] <- "parser_failed"
    row$adjudication_basis[
      parser_failed
    ] <- "The parser did not produce a glycan."

    ref_match <- row$reference_comparison_status == "equivalent"
    ref_diff <- row$reference_comparison_status == "semantic_difference"
    row$adjudication_status[ref_match] <- "validated_by_accession_reference"
    row$adjudication_basis[ref_match] <- paste0(
      "The glyparse canonical glycan equals the independently normalized ",
      "accession-matched IUPAC-condensed reference."
    )

    for (index in which(!parser_failed & !ref_match)) {
      available <- nzchar(external_matrix[index, ])
      values <- external_matrix[index, available]
      scopes <- scope_matrix[index, available]
      direct <- grepl("^direct", scopes)
      direct_parser_hits <- sum(
        direct & values == row$parsed_iupac_condensed[[index]]
      )
      reference <- row$reference_normalized_iupac_condensed[[index]]
      direct_reference_hits <- if (nzchar(reference)) {
        sum(direct & values == reference)
      } else {
        0L
      }
      direct_values <- values[direct]
      direct_counts <- sort(table(direct_values), decreasing = TRUE)
      direct_consensus <- length(direct_counts) > 0L &&
        direct_counts[[1L]] >= 2L
      direct_consensus_value <- if (direct_consensus) {
        names(direct_counts)[[1L]]
      } else {
        ""
      }

      if (ref_diff[[index]] && direct_reference_hits > 0L) {
        row$adjudication_status[[
          index
        ]] <- "glyparse_reference_difference_supported"
        row$adjudication_basis[[index]] <- paste0(
          "glyparse differs from the accession reference, and ",
          direct_reference_hits,
          " direct normalized external result(s) reproduce the reference glycan."
        )
      } else if (
        ref_diff[[index]] &&
          direct_parser_hits > 0L &&
          direct_reference_hits == 0L
      ) {
        row$adjudication_status[[index]] <- "reference_difference_supported"
        row$adjudication_basis[[index]] <- paste0(
          "glyparse differs from the accession reference, but ",
          direct_parser_hits,
          " direct normalized external result(s) reproduce the glyparse glycan."
        )
      } else if (ref_diff[[index]]) {
        row$adjudication_status[[index]] <- "unresolved_reference_difference"
        row$adjudication_basis[[index]] <- paste0(
          "glyparse differs from the accession reference; available external ",
          "results do not resolve which representation is correct."
        )
      } else if (
        direct_consensus &&
          !identical(
            direct_consensus_value,
            row$parsed_iupac_condensed[[index]]
          )
      ) {
        row$adjudication_status[[
          index
        ]] <- "unresolved_external_consensus_difference"
        row$adjudication_basis[[index]] <- paste0(
          "At least two direct external converters agree with each other and ",
          "differ from glyparse, but shared external information loss cannot ",
          "be excluded without an accession reference."
        )
      } else if (direct_parser_hits > 0L) {
        row$adjudication_status[[index]] <- "supported_by_external_result"
        row$adjudication_basis[[index]] <- paste0(
          direct_parser_hits,
          " direct normalized external result(s) reproduce the glyparse glycan."
        )
      } else if (length(values) > 0L) {
        row$adjudication_status[[index]] <- "unresolved_external_difference"
        row$adjudication_basis[[index]] <- paste0(
          "External normalized result(s) are available but differ from glyparse ",
          "without an independent consensus sufficient for attribution."
        )
      }
    }
  }

  failures <- row[
    row$parse_status == "failed",
    c(
      "format",
      "format_label",
      "parser",
      "row_index",
      "glytoucan_ac",
      "source_sequence",
      "parse_error"
    )
  ]
  if (nrow(failures) > 0L) {
    failure_results[[format]] <- cbind(
      failures,
      diagnose_failures(format, failures)
    )
  }

  row_results[[format]] <- row
  format_summaries[[format]] <- data.frame(
    format = format,
    format_label = unique(parser$format_label),
    rows = n,
    parsed = sum(parser$parse_status == "parsed"),
    failed = sum(parser$parse_status == "failed"),
    parse_success_rate = mean(parser$parse_status == "parsed"),
    reference_available = sum(parser$reference_validation_applicable),
    reference_equivalent = sum(
      parser$reference_comparison_status == "equivalent"
    ),
    reference_semantic_difference = sum(
      parser$reference_comparison_status == "semantic_difference"
    ),
    validated_or_supported = sum(
      row$adjudication_status %in%
        c(
          "validated_by_accession_reference",
          "supported_by_external_result",
          "reference_difference_supported"
        )
    ),
    glyparse_outlier = sum(
      row$adjudication_status %in%
        c("glyparse_reference_difference_supported")
    ),
    unresolved = sum(grepl("^unresolved", row$adjudication_status)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

rbind_fill <- function(tables) {
  fields <- Reduce(union, lapply(tables, names))
  tables <- lapply(tables, function(table) {
    missing <- setdiff(fields, names(table))
    for (field in missing) {
      table[[field]] <- ""
    }
    table[, fields, drop = FALSE]
  })
  do.call(rbind, tables)
}

row_results <- rbind_fill(row_results)
failure_results <- do.call(rbind, failure_results)
difference_results <- if (length(difference_results) == 0L) {
  data.frame()
} else {
  do.call(rbind, difference_results)
}
format_summaries <- do.call(rbind, format_summaries)
tool_summaries <- do.call(rbind, tool_summaries)
rownames(row_results) <- NULL
rownames(failure_results) <- NULL
rownames(difference_results) <- NULL
rownames(format_summaries) <- NULL
rownames(tool_summaries) <- NULL

benchmark_write_csv(
  row_results,
  file.path(benchmark_output_dir, "validation_rows.csv.gz")
)
benchmark_write_csv(
  failure_results,
  file.path(benchmark_output_dir, "parse_failures.csv.gz")
)
benchmark_write_csv(
  difference_results,
  file.path(benchmark_output_dir, "semantic_disagreements.csv.gz")
)
benchmark_write_csv(
  format_summaries,
  file.path(benchmark_output_dir, "format_summary.csv")
)
benchmark_write_csv(
  tool_summaries,
  file.path(benchmark_output_dir, "external_tool_summary.csv")
)

benchmark_repeats <- benchmark_read_csv(file.path(
  benchmark_cache_dir,
  "glyparse-benchmark-repeats.csv"
))
benchmark_summary <- do.call(
  rbind,
  lapply(
    split(benchmark_repeats, benchmark_repeats$format),
    function(group) {
      data.frame(
        format = group$format[[1L]],
        parser = group$parser[[1L]],
        mono_type = group$mono_type[[1L]],
        sample_rows = group$sample_rows[[1L]],
        repeats = nrow(group),
        median_elapsed_seconds = stats::median(group$elapsed_seconds),
        median_rows_per_second = stats::median(group$rows_per_second),
        min_rows_per_second = min(group$rows_per_second),
        max_rows_per_second = max(group$rows_per_second),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  )
)
rownames(benchmark_summary) <- NULL
benchmark_write_csv(
  benchmark_summary,
  file.path(benchmark_output_dir, "parser_benchmark.csv")
)

failure_reason_summary <- aggregate(
  list(rows = rep(1L, nrow(failure_results))),
  by = failure_results[, c(
    "format",
    "failure_class",
    "failure_stage",
    "failure_reason"
  )],
  FUN = sum
)
failure_reason_summary <- failure_reason_summary[
  order(failure_reason_summary$format, -failure_reason_summary$rows),
]
benchmark_write_csv(
  failure_reason_summary,
  file.path(benchmark_output_dir, "failure_reason_summary.csv")
)

stopifnot(
  nrow(row_results) == sum(format_summaries$rows),
  nrow(failure_results) == sum(format_summaries$failed),
  all(nzchar(failure_results$failure_reason)),
  all(nzchar(failure_results$failure_class)),
  all(
    !grepl("^glyparse_.*outlier$", row_results$adjudication_status) |
      row_results$parse_status == "parsed"
  )
)

print(format_summaries)
print(tool_summaries)
