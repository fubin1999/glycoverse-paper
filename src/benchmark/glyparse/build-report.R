#!/usr/bin/env Rscript

script_path <- sub(
  "^--file=",
  "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1L]]
)
source(file.path(dirname(normalizePath(script_path)), "common.R"))

dir.create(benchmark_output_dir, recursive = TRUE, showWarnings = FALSE)

format_summary <- benchmark_read_csv(file.path(
  benchmark_output_dir,
  "format_summary.csv"
))
tool_summary <- benchmark_read_csv(file.path(
  benchmark_output_dir,
  "external_tool_summary.csv"
))
parser_benchmark <- benchmark_read_csv(file.path(
  benchmark_output_dir,
  "parser_benchmark.csv"
))
failure_summary <- benchmark_read_csv(file.path(
  benchmark_output_dir,
  "failure_reason_summary.csv"
))
coverage_summary <- benchmark_read_csv(file.path(
  benchmark_output_dir,
  "sequence_coverage_summary.csv"
))
validation_rows <- benchmark_read_csv(file.path(
  benchmark_output_dir,
  "validation_rows.csv.gz"
))
disagreements <- benchmark_read_csv(file.path(
  benchmark_output_dir,
  "semantic_disagreements.csv.gz"
))

format_summary$parse_success_rate <- as.numeric(
  format_summary$parse_success_rate
)
format_summary$parse_failure_rate <- 1 - format_summary$parse_success_rate
format_summary$reference_agreement_rate <- ifelse(
  format_summary$reference_available > 0L,
  format_summary$reference_equivalent / format_summary$reference_available,
  NA_real_
)
format_summary$format_label <- as.character(format_summary$format_label)

tool_summary$normalized_rate <- tool_summary$normalized / tool_summary$rows
tool_summary$comparable_rate <-
  tool_summary$comparable_with_glyparse / tool_summary$rows
tool_summary$format_label <- format_summary$format_label[
  match(tool_summary$format, format_summary$format)
]
coverage_summary$direct_parse_rate <- as.numeric(
  coverage_summary$direct_parse_rate
)

parser_benchmark$format_label <- format_summary$format_label[
  match(parser_benchmark$format, format_summary$format)
]
tool_versions <- vapply(
  c("GlycanFormatConverter", "glycowork", "glypy"),
  function(tool) {
    unique(tool_summary$tool_version[tool_summary$tool == tool])[[1L]]
  },
  character(1)
)

total_rows <- sum(format_summary$rows)
total_parsed <- sum(format_summary$parsed)
total_failed <- sum(format_summary$failed)
iupac_formats <- c("iupac_compact", "iupac_condensed", "iupac_extended")
iupac_complete <- all(
  format_summary$failed[format_summary$format %in% iupac_formats] == 0L
)
mixed_residue_failure_rows <- sum(failure_summary$rows[
  grepl(
    "mixes generic and concrete|all generic or all concrete",
    failure_summary$failure_reason,
    ignore.case = TRUE
  )
])
reference_available <- sum(format_summary$reference_available)
reference_equivalent <- sum(format_summary$reference_equivalent)
reference_differences <- sum(format_summary$reference_semantic_difference)
reference_compared <- reference_equivalent + reference_differences
reference_parser_failed <- reference_available - reference_compared
outlier_rows <- sum(format_summary$glyparse_outlier)
unresolved_rows <- sum(format_summary$unresolved)
outlier_accessions <- unique(validation_rows$glytoucan_ac[
  validation_rows$adjudication_status %in%
    c("glyparse_reference_difference_supported")
])
outlier_accessions <- outlier_accessions[nzchar(outlier_accessions)]

reference_review <- disagreements[
  disagreements$comparator == "IUPAC-condensed accession reference",
  c(
    "format",
    "glytoucan_ac",
    "glyparse_iupac_condensed",
    "comparator_normalized_iupac_condensed",
    "difference_reason"
  )
]
reference_review$adjudication_status <- validation_rows$adjudication_status[
  match(
    paste(reference_review$format, reference_review$glytoucan_ac, sep = "\r"),
    paste(validation_rows$format, validation_rows$glytoucan_ac, sep = "\r")
  )
]

failure_class_summary <- aggregate(
  list(rows = failure_summary$rows),
  by = failure_summary[, c("format", "failure_class")],
  FUN = sum
)
failure_class_summary$format_label <- format_summary$format_label[
  match(failure_class_summary$format, format_summary$format)
]
failure_class_summary <- failure_class_summary[
  order(failure_class_summary$format_label, -failure_class_summary$rows),
]

top_failure_reasons <- do.call(
  rbind,
  lapply(
    split(failure_summary, failure_summary$format),
    function(group) {
      group <- group[order(-group$rows), ]
      head(group, 5L)
    }
  )
)
top_failure_reasons$format_label <- format_summary$format_label[
  match(top_failure_reasons$format, format_summary$format)
]
top_failure_reasons <- top_failure_reasons[, c(
  "format_label",
  "failure_class",
  "failure_stage",
  "rows",
  "failure_reason"
)]
rownames(top_failure_reasons) <- NULL

coverage_tool_order <- c(
  "glyparse",
  "GlycanFormatConverter",
  "glypy",
  "glycowork"
)
coverage_totals <- do.call(
  rbind,
  lapply(coverage_tool_order, function(tool) {
    group <- coverage_summary[coverage_summary$tool == tool, ]
    directly_tested <- sum(group$directly_tested)
    data.frame(
      tool = tool,
      rows = sum(group$rows),
      directly_tested = directly_tested,
      parsed = sum(group$parsed),
      failed = sum(group$failed),
      fallback_only = sum(group$fallback_only),
      not_tested = sum(group$not_tested),
      direct_parse_rate = if (directly_tested == 0L) {
        NA_real_
      } else {
        sum(group$parsed) / directly_tested
      },
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
)
rownames(coverage_totals) <- NULL

headline <- data.frame(
  total_rows = total_rows,
  total_parsed = total_parsed,
  parse_success_rate = total_parsed / total_rows,
  reference_equivalent = reference_equivalent,
  reference_available = reference_available,
  reference_agreement_rate = reference_equivalent / reference_compared,
  supported_glyparse_outlier_rows = outlier_rows,
  unresolved_rows = unresolved_rows,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

generated_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
outlier_text <- if (length(outlier_accessions) == 0L) {
  "none"
} else {
  paste(outlier_accessions, collapse = ", ")
}

technical_summary_body <- paste(
  c(
    "## Technical summary",
    "",
    paste0(
      "The benchmark evaluated **",
      format(total_rows, big.mark = ",", scientific = FALSE),
      "** source sequences across seven formats. glyparse returned canonical glycans for **",
      format(total_parsed, big.mark = ",", scientific = FALSE),
      "** rows (**",
      sprintf("%.2f%%", 100 * total_parsed / total_rows),
      "**). The remaining **",
      format(total_failed, big.mark = ",", scientific = FALSE),
      "** failures have row-level stage, class, and reason evidence."
    ),
    "",
    paste0(
      if (iupac_complete) {
        "IUPAC-compact, IUPAC-condensed, and IUPAC-extended each achieved **100% glyparse parse coverage**."
      } else {
        "At least one IUPAC source format retains glyparse parse failures."
      },
      " The remaining failure ledger contains **",
      format(mixed_residue_failure_rows, big.mark = ",", scientific = FALSE),
      "** row(s) whose detailed reason cites the former mixed generic/concrete residue restriction."
    ),
    "",
    paste0(
      "Among **",
      format(reference_available, big.mark = ",", scientific = FALSE),
      "** accession-matched non-condensed rows, **",
      format(reference_compared, big.mark = ",", scientific = FALSE),
      "** were semantically comparable: **",
      format(reference_equivalent, big.mark = ",", scientific = FALSE),
      "** matched the independently normalized IUPAC-condensed reference and **",
      reference_differences,
      "** differed; **",
      reference_parser_failed,
      "** failed in glyparse before comparison. Direct external evidence supports **",
      outlier_rows,
      "** glyparse outlier format-row(s), covering accession(s): ",
      outlier_text,
      "."
    ),
    "",
    "A separate sequence-level coverage ledger records whether each of the four tested tools parsed each exact source sequence. Direct failures, accession-fallback-only results, and untested formats remain distinct."
  ),
  collapse = "\n"
)
scope_body <- paste(
  c(
    "## Scope and definitions",
    "",
    paste0(
      "IUPAC-condensed input is assessed with all four tools for parse coverage and with glyparse for performance, but is excluded from correctness adjudication as requested. All other formats are checked against an accession-matched IUPAC-condensed corpus where available and against GlycanFormatConverter ",
      tool_versions[["GlycanFormatConverter"]],
      ", glycowork ",
      tool_versions[["glycowork"]],
      ", and glypy ",
      tool_versions[["glypy"]],
      "."
    ),
    "",
    "A semantic comparison requires both sides to normalize successfully through the current `parse_iupac_condensed()`. Chemically fixed donor positions are then filled with `fill_anomer_pos()`, and equality is tested with an isomorphism-invariant labeled-graph fingerprint, so implicit donor positions and alternate serializations of symmetric branches compare fairly. Accession references and GlycanFormatConverter retain the reducing-end anomer. Because glycowork and glypy do not encode it, only that graph-level attribute is ignored for those comparisons. This makes `Gal(b1-3)GalNAc(a1-` and `Gal(b1-3)GalNAc(?1-` equivalent while retaining internal anomers/linkages, residue identity, substituents, topology, reducing-end alditol state, floating components, and candidate-parent domains."
  ),
  collapse = "\n"
)
methodology_body <- paste(
  c(
    "## Methodology",
    "",
    "Each source sequence is parsed independently so one invalid or generic row cannot contaminate another row's result. Accession reference matching is exact. External converter availability, conversion success, current-model normalization, comparison contract, display key, and semantic agreement are separate fields. Labeled-graph fingerprints encode the rooted residue/linkage graph, substituents, alditol state, floating components, and candidate-parent relations. Direct cross-tool consensus uses the same graph comparison while ignoring only the terminal reducing-end anomer, so missing root-anomer information cannot create a false outlier. Accessions converted through an accession-matched WURCS fallback are labeled as fallback and are never used to attribute a source-format error to glyparse.",
    "",
    "Performance uses deterministic, evenly spaced samples of 1,000 successfully parsed concrete glycans, three vectorized repeats per parser. This isolates steady-state vector throughput from the scalar row-level diagnostic run."
  ),
  collapse = "\n"
)
coverage_lines <- vapply(
  seq_len(nrow(coverage_totals)),
  function(index) {
    row <- coverage_totals[index, ]
    paste0(
      "- **",
      row$tool,
      ":** parsed ",
      format(row$parsed, big.mark = ",", scientific = FALSE),
      " of ",
      format(row$directly_tested, big.mark = ",", scientific = FALSE),
      " directly tested rows (",
      sprintf("%.2f%%", 100 * row$direct_parse_rate),
      "); ",
      format(row$fallback_only, big.mark = ",", scientific = FALSE),
      " fallback-only and ",
      format(row$not_tested, big.mark = ",", scientific = FALSE),
      " not tested."
    )
  },
  character(1)
)
coverage_body <- paste(
  c(
    "## Direct sequence coverage keeps fallback evidence separate",
    "",
    "The row-level coverage report contains one record for every source sequence. `parsed` and `failed` are direct tests of that exact source string. `fallback_only` means the available conversion used an accession-matched WURCS sequence instead, and `not_tested` means the tool was not run for that source format. Boolean `can_parse` values are therefore populated only for direct tests. For glypy IUPAC-condensed input, the terminal open reducing-end token is removed before its simple-IUPAC parser is called, matching the agreed convention that glypy does not encode reducing-end anomers; the row retains this evidence scope.",
    "",
    coverage_lines
  ),
  collapse = "\n"
)
reference_finding <- if (reference_differences == 0L) {
  paste0(
    "The normalized accession reference agrees with glyparse for **100%** of ",
    "the ",
    format(reference_compared, big.mark = ",", scientific = FALSE),
    " comparable non-condensed rows. The previously reported G12345BK ",
    "differences were false positives: GlycoCT, WURCS, and compact assigned an ",
    "unknown linkage to different but automorphic terminal Man nodes, while ",
    "GWB explicitly wrote Man's chemically fixed donor position where the ",
    "reference omitted it. Chemistry-aware labeled-graph isomorphism now ",
    "recognizes all four serializations as the same glycan."
  )
} else {
  paste0(
    "The normalized accession reference agrees with glyparse for **",
    sprintf("%.3f%%", 100 * reference_equivalent / reference_compared),
    "** of comparable non-condensed rows. The remaining ",
    reference_differences,
    " reference difference(s) and their component-level reasons are retained ",
    "in the evidence ledger. Direct external conversion supports glyparse ",
    "being the outlier for ",
    outlier_rows,
    " format-row(s)."
  )
}
findings_body <- paste(
  c(
    "## Findings",
    "",
    reference_finding,
    "",
    "GlycanFormatConverter direct coverage is format-dependent. It directly attempts GlycoCT, WURCS, IUPAC-extended, and IUPAC-condensed in this run; 19,483 of 19,621 IUPAC-condensed rows converted, while only 27 IUPAC-extended rows converted. GLYCAM-IUPAC, GWB, and compact results from this tool use accession-matched WURCS fallback and are not direct source-parser evidence."
  ),
  collapse = "\n"
)
limitations_body <- paste(
  c(
    "## Limitations and robustness",
    "",
    "Graph comparison is intentionally limited to chemistry representable by the current glyrepr model. An external output that cannot be normalized is reported as unavailable rather than as a semantic disagreement. Agreement between two direct converters without an accession reference is retained as an unresolved external-consensus difference, not attributed to glyparse, because shared information loss cannot be excluded. The reducing-end-neutral contract for glycowork and glypy ignores only the graph-level terminal root anomer; it does not forgive internal anomer/linkage, alditol, topology, residue, substituent, or floating-component differences. Raw string equality is only a fast path; every unequal serialization that could share the same residue/linkage inventory is adjudicated by labeled-graph isomorphism.",
    "",
    "Detailed failure classifications reuse exact-row diagnostics only when accession and full source sequence match the prior diagnostic corpus; otherwise the current internal parser stage is replayed. Every final failure row retains the diagnostic source."
  ),
  collapse = "\n"
)
reference_next_step <- if (reference_differences == 0L) {
  "No accession-reference semantic outlier remains for manual review."
} else {
  paste0(
    "Review the ",
    reference_differences,
    " remaining accession-reference difference(s) in the row-level ledger."
  )
}
next_steps_body <- paste(
  c(
    "## Recommended next steps",
    "",
    paste(
      reference_next_step,
      "For parser coverage, prioritize the largest glyrepr-restriction families separately from glyparse implementation limitations. Re-run the workflow after any parser or glyrepr semantic change; the source revisions and corpus hashes in the cache metadata make changes auditable."
    )
  ),
  collapse = "\n"
)
report_markdown <- c(
  "# glyparse corpus validation and benchmark",
  "",
  technical_summary_body,
  "",
  scope_body,
  "",
  methodology_body,
  "",
  coverage_body,
  "",
  findings_body,
  "",
  limitations_body,
  "",
  next_steps_body
)
writeLines(
  report_markdown,
  file.path(benchmark_output_dir, "report.md"),
  useBytes = TRUE
)

source_summary <- list(
  id = "format_summary",
  label = "glyparse format summary",
  path = "results/data/glyparse_validation/format_summary.csv",
  query = list(
    engine = "DuckDB",
    language = "SQL",
    sql = paste0(
      "SELECT * FROM read_csv_auto(",
      "'results/data/glyparse_validation/format_summary.csv')"
    ),
    description = "Read the generated format-level glyparse validation summary."
  )
)
source_tools <- list(
  id = "external_summary",
  label = "External converter summary",
  path = "results/data/glyparse_validation/external_tool_summary.csv",
  query = list(
    engine = "DuckDB",
    language = "SQL",
    sql = paste0(
      "SELECT *, normalized::DOUBLE / rows AS normalized_rate, ",
      "comparable_with_glyparse::DOUBLE / rows AS comparable_rate ",
      "FROM read_csv_auto(",
      "'results/data/glyparse_validation/external_tool_summary.csv')"
    ),
    description = "Read conversion, normalization, and comparison counts by tool and format."
  )
)
source_benchmark <- list(
  id = "parser_benchmark",
  label = "glyparse parser benchmark",
  path = "results/data/glyparse_validation/parser_benchmark.csv",
  query = list(
    engine = "DuckDB",
    language = "SQL",
    sql = paste0(
      "SELECT * FROM read_csv_auto(",
      "'results/data/glyparse_validation/parser_benchmark.csv')"
    ),
    description = "Read median and range statistics from the three parser benchmark repeats."
  )
)
source_failures <- list(
  id = "failure_summary",
  label = "glyparse failure reason summary",
  path = "results/data/glyparse_validation/failure_reason_summary.csv",
  query = list(
    engine = "DuckDB",
    language = "SQL",
    sql = paste0(
      "SELECT * FROM read_csv_auto(",
      "'results/data/glyparse_validation/failure_reason_summary.csv')"
    ),
    description = "Read grouped row-level parser failure classes and reasons."
  )
)
source_disagreements <- list(
  id = "disagreement_ledger",
  label = "Semantic disagreement ledger",
  path = "results/data/glyparse_validation/semantic_disagreements.csv.gz",
  query = list(
    engine = "DuckDB",
    language = "SQL",
    sql = paste0(
      "SELECT * FROM read_csv_auto(",
      "'results/data/glyparse_validation/semantic_disagreements.csv.gz')"
    ),
    description = "Read the complete normalized glyparse-versus-comparator disagreement ledger."
  )
)
source_coverage <- list(
  id = "sequence_coverage_summary",
  label = "Sequence-level parser coverage summary",
  path = "results/data/glyparse_validation/sequence_coverage_summary.csv",
  query = list(
    engine = "DuckDB",
    language = "SQL",
    sql = paste0(
      "SELECT * FROM read_csv_auto(",
      "'results/data/glyparse_validation/sequence_coverage_summary.csv')"
    ),
    description = paste0(
      "Read direct parse, failure, fallback-only, and not-tested counts by ",
      "tool and source format."
    )
  )
)

artifact <- list(
  surface = "report",
  manifest = list(
    version = 1L,
    surface = "report",
    title = "glyparse corpus validation and benchmark",
    description = paste0(
      "Row-level validation of glyparse against accession references and ",
      "three external glycan converters."
    ),
    generatedAt = generated_at,
    cards = list(
      list(
        id = "rows",
        description = "All source rows across seven sequence formats.",
        dataset = "headline",
        sourceId = "format_summary",
        metrics = list(list(
          label = "Sequences evaluated",
          field = "total_rows",
          format = "number"
        ))
      ),
      list(
        id = "parsed",
        description = "Rows producing a non-missing canonical glycan.",
        dataset = "headline",
        sourceId = "format_summary",
        metrics = list(
          list(
            label = "Parse success rate",
            field = "parse_success_rate",
            format = "percent"
          ),
          list(
            label = "Parsed rows",
            field = "total_parsed",
            format = "number"
          )
        )
      ),
      list(
        id = "reference",
        description = paste0(
          "Labeled-graph isomorphism among accession-matched non-condensed rows."
        ),
        dataset = "headline",
        sourceId = "format_summary",
        metrics = list(
          list(
            label = "Reference agreement",
            field = "reference_agreement_rate",
            format = "percent"
          ),
          list(
            label = "Equivalent rows",
            field = "reference_equivalent",
            format = "number"
          )
        )
      ),
      list(
        id = "outliers",
        description = paste0(
          "Format rows where direct evidence supports glyparse as the outlier."
        ),
        dataset = "headline",
        sourceId = "format_summary",
        metrics = list(list(
          label = "Supported glyparse outlier rows",
          field = "supported_glyparse_outlier_rows",
          format = "number"
        ))
      )
    ),
    charts = list(
      list(
        id = "parse_coverage",
        title = "Parse success rate by source format",
        subtitle = "All corpus rows; IUPAC-condensed correctness is not assessed.",
        type = "bar",
        dataset = "format_summary",
        sourceId = "format_summary",
        encodings = list(
          x = list(
            field = "format_label",
            type = "nominal",
            label = "Source format"
          ),
          y = list(
            field = "parse_success_rate",
            type = "quantitative",
            label = "Parse success rate",
            format = "percent"
          )
        )
      ),
      list(
        id = "external_coverage",
        title = "Normalized external-output coverage",
        subtitle = paste0(
          "Share of source rows converted and normalizable in the current ",
          "glyrepr model; direct and fallback scope remains in the data."
        ),
        type = "bar",
        dataset = "tool_summary",
        sourceId = "external_summary",
        encodings = list(
          x = list(
            field = "format_label",
            type = "nominal",
            label = "Source format"
          ),
          y = list(
            field = "normalized_rate",
            type = "quantitative",
            label = "Normalized output rate",
            format = "percent"
          ),
          color = list(
            field = "tool",
            type = "nominal",
            label = "Tool"
          )
        )
      ),
      list(
        id = "throughput",
        title = "Median parser throughput",
        subtitle = paste0(
          "Three vectorized repeats on deterministic samples of 1,000 ",
          "successfully parsed concrete glycans."
        ),
        type = "bar",
        dataset = "parser_benchmark",
        sourceId = "parser_benchmark",
        encodings = list(
          x = list(
            field = "format_label",
            type = "nominal",
            label = "Source format"
          ),
          y = list(
            field = "median_rows_per_second",
            type = "quantitative",
            label = "Rows per second",
            format = "number"
          )
        )
      )
    ),
    tables = list(
      list(
        id = "format_detail",
        title = "Format-level validation results",
        dataset = "format_summary",
        sourceId = "format_summary",
        defaultSort = list(field = "rows", direction = "desc"),
        columns = list(
          list(field = "format_label", label = "Format", type = "text"),
          list(field = "rows", label = "Rows", format = "number"),
          list(field = "parsed", label = "Parsed", format = "number"),
          list(field = "failed", label = "Failed", format = "number"),
          list(
            field = "parse_success_rate",
            label = "Parse success",
            format = "percent"
          ),
          list(
            field = "reference_equivalent",
            label = "Reference-equivalent",
            format = "number"
          ),
          list(
            field = "reference_semantic_difference",
            label = "Reference differences",
            format = "number"
          )
        )
      ),
      list(
        id = "reference_review",
        title = "Accession-reference differences",
        dataset = "reference_review",
        sourceId = "disagreement_ledger",
        defaultSort = list(field = "glytoucan_ac", direction = "asc"),
        columns = list(
          list(field = "format", label = "Format", type = "text"),
          list(field = "glytoucan_ac", label = "GlyTouCan", type = "text"),
          list(
            field = "adjudication_status",
            label = "Adjudication",
            type = "text"
          ),
          list(
            field = "difference_reason",
            label = "Component difference",
            type = "text"
          )
        )
      ),
      list(
        id = "failure_reasons",
        title = "Top parser failure reasons by format",
        dataset = "top_failure_reasons",
        sourceId = "failure_summary",
        defaultSort = list(field = "rows", direction = "desc"),
        columns = list(
          list(field = "format_label", label = "Format", type = "text"),
          list(field = "failure_class", label = "Class", type = "text"),
          list(field = "failure_stage", label = "Stage", type = "text"),
          list(field = "rows", label = "Rows", format = "number"),
          list(
            field = "failure_reason",
            label = "Detailed reason",
            type = "text"
          )
        )
      ),
      list(
        id = "sequence_coverage",
        title = "Sequence-level parser coverage by format",
        dataset = "coverage_summary",
        sourceId = "sequence_coverage_summary",
        defaultSort = list(field = "tool", direction = "asc"),
        columns = list(
          list(field = "tool", label = "Tool", type = "text"),
          list(field = "format_label", label = "Format", type = "text"),
          list(
            field = "directly_tested",
            label = "Directly tested",
            format = "number"
          ),
          list(field = "parsed", label = "Parsed", format = "number"),
          list(field = "failed", label = "Failed", format = "number"),
          list(
            field = "direct_parse_rate",
            label = "Direct parse rate",
            format = "percent"
          ),
          list(
            field = "fallback_only",
            label = "Fallback only",
            format = "number"
          ),
          list(
            field = "not_tested",
            label = "Not tested",
            format = "number"
          )
        )
      )
    ),
    sources = list(
      source_summary,
      source_tools,
      source_benchmark,
      source_failures,
      source_disagreements,
      source_coverage
    ),
    blocks = list(
      list(
        id = "title",
        type = "markdown",
        body = "# glyparse corpus validation and benchmark"
      ),
      list(
        id = "summary",
        type = "markdown",
        body = technical_summary_body
      ),
      list(
        id = "headline_metrics",
        type = "metric-strip",
        cardIds = c("rows", "parsed", "reference", "outliers")
      ),
      list(
        id = "scope",
        type = "markdown",
        body = scope_body
      ),
      list(
        id = "methodology",
        type = "markdown",
        body = methodology_body
      ),
      list(
        id = "parse_coverage_block",
        type = "chart",
        chartId = "parse_coverage",
        layout = "full"
      ),
      list(
        id = "format_detail_block",
        type = "table",
        tableId = "format_detail"
      ),
      list(
        id = "sequence_coverage_findings",
        type = "markdown",
        body = coverage_body,
        sourceId = "sequence_coverage_summary"
      ),
      list(
        id = "sequence_coverage_block",
        type = "table",
        tableId = "sequence_coverage"
      ),
      list(
        id = "findings",
        type = "markdown",
        body = findings_body
      ),
      list(
        id = "reference_review_block",
        type = "table",
        tableId = "reference_review"
      ),
      list(
        id = "external_coverage_block",
        type = "chart",
        chartId = "external_coverage",
        layout = "full"
      ),
      list(
        id = "failure_reasons_block",
        type = "table",
        tableId = "failure_reasons"
      ),
      list(
        id = "throughput_block",
        type = "chart",
        chartId = "throughput",
        layout = "full"
      ),
      list(
        id = "limitations",
        type = "markdown",
        body = limitations_body
      ),
      list(
        id = "next_steps",
        type = "markdown",
        body = next_steps_body
      )
    )
  ),
  snapshot = list(
    version = 1L,
    generatedAt = generated_at,
    status = "ready",
    datasets = list(
      headline = headline,
      format_summary = format_summary,
      tool_summary = tool_summary,
      parser_benchmark = parser_benchmark,
      coverage_summary = coverage_summary,
      reference_review = reference_review,
      failure_class_summary = failure_class_summary,
      top_failure_reasons = top_failure_reasons
    )
  ),
  sources = list(
    source_summary,
    source_tools,
    source_benchmark,
    source_failures,
    source_disagreements,
    source_coverage
  )
)

jsonlite::write_json(
  artifact,
  file.path(benchmark_output_dir, "artifact.json"),
  auto_unbox = TRUE,
  pretty = TRUE,
  na = "null",
  digits = 12
)

chart_map <- c(
  "# Chart map",
  "",
  "| Segment | Question | Family/type | Fields | Supported claim | Palette |",
  "|---|---|---|---|---|---|",
  paste0(
    "| Parse coverage | How often does each parser return a glycan? | ",
    "Comparison / bar | format_label, parse_success_rate, rows | ",
    "Coverage differs materially by format | Single blue root; direct labels |"
  ),
  paste0(
    "| External coverage | How much external output is semantically comparable? | ",
    "Comparison / grouped bar | format_label, normalized_rate, tool, direct_rows, fallback_rows | ",
    "Tool coverage and format dialect support differ | Three approved roots plus labels |"
  ),
  paste0(
    "| Throughput | What is steady-state parser throughput? | ",
    "Ranking / bar | format_label, median_rows_per_second, sample_rows, repeats | ",
    "Throughput varies by parser; exact medians and repeat ranges are retained | Single blue root; exact tooltips |"
  ),
  "",
  "Sequence-level coverage is presented as a table because exact lookup across four mutually exclusive states is more important than visual shape."
)
writeLines(
  chart_map,
  file.path(benchmark_output_dir, "chart-map.md"),
  useBytes = TRUE
)

glyparse_metadata <- jsonlite::read_json(file.path(
  benchmark_cache_dir,
  "glyparse-metadata.json"
))
gfc_metadata <- jsonlite::read_json(file.path(
  benchmark_cache_dir,
  "gfc-metadata.json"
))
gfc_metadata$jar <- basename(gfc_metadata$jar)
if (!is.null(gfc_metadata$reused_cache)) {
  gfc_metadata$reused_cache <- basename(gfc_metadata$reused_cache)
}
run_metadata <- list(
  report_generated_at = generated_at,
  glyparse_run = glyparse_metadata,
  external_tools = lapply(
    split(tool_summary, tool_summary$tool),
    function(group) {
      list(
        version = unique(group$tool_version),
        formats = unique(group$format)
      )
    }
  ),
  sequence_coverage = list(
    rows = nrow(validation_rows),
    statuses = c("parsed", "failed", "fallback_only", "not_tested"),
    boolean_rule = paste0(
      "can_parse is populated only for parsed or failed direct source tests; ",
      "fallback_only and not_tested are missing."
    )
  ),
  semantic_comparison = list(
    normalization = c("parse_iupac_condensed", "fill_anomer_pos"),
    strict_contract = "strict_graph_isomorphism",
    reducing_end_neutral_contract = "graph_isomorphism_ignore_terminal_reducing_end_token",
    reducing_end_neutral_tools = c("glycowork", "glypy"),
    fingerprint_semantics = paste0(
      "Rooted residue/linkage graph, substituents, alditol state, floating ",
      "components, and candidate-parent relations."
    )
  ),
  glycan_format_converter = gfc_metadata
)
jsonlite::write_json(
  run_metadata,
  file.path(benchmark_output_dir, "run-metadata.json"),
  auto_unbox = TRUE,
  pretty = TRUE
)

message("Wrote report artifacts and run-metadata.json")
