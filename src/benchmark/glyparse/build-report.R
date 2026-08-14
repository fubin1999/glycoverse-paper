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

parser_benchmark$format_label <- format_summary$format_label[
  match(parser_benchmark$format, format_summary$format)
]

total_rows <- sum(format_summary$rows)
total_parsed <- sum(format_summary$parsed)
total_failed <- sum(format_summary$failed)
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
    )
  ),
  collapse = "\n"
)
scope_body <- paste(
  c(
    "## Scope and definitions",
    "",
    "IUPAC-condensed input is assessed only for parse coverage and performance, as requested. All other formats are checked against an accession-matched IUPAC-condensed corpus where available and against GlycanFormatConverter 2.10.3, glycowork 1.8.0, and glypy 1.0.17.",
    "",
    "A semantic comparison requires both sides to normalize successfully through the current `parse_iupac_condensed()` and then have identical glyrepr canonical serialization. This removes branch and node-order differences while preserving residue identity, substituents, topology, linkages, reducing-end anomer and alditol state, floating components, and candidate-parent domains."
  ),
  collapse = "\n"
)
methodology_body <- paste(
  c(
    "## Methodology",
    "",
    "Each source sequence is parsed independently so one invalid or generic row cannot contaminate another row's result. Accession reference matching is exact. External converter availability, conversion success, current-model normalization, and semantic agreement are separate fields. Accessions converted through an accession-matched WURCS fallback are labeled as fallback and are never used to attribute a source-format error to glyparse.",
    "",
    "Performance uses deterministic, evenly spaced samples of 1,000 successfully parsed concrete glycans, three vectorized repeats per parser. This isolates steady-state vector throughput from the scalar row-level diagnostic run."
  ),
  collapse = "\n"
)
findings_body <- paste(
  c(
    "## Findings",
    "",
    paste0(
      "The normalized accession reference agrees with glyparse for **",
      sprintf("%.3f%%", 100 * reference_equivalent / reference_compared),
      "** of comparable non-condensed rows. The five reference differences represent two accessions and concern anomeric-carbon specificity on mannose residues. Direct external conversion supports glyparse being the outlier for ",
      outlier_rows,
      " format-row(s); the other reference differences remain unresolved or support the reference/source distinction described in the evidence table."
    ),
    "",
    "GlycanFormatConverter direct coverage is format-dependent. It directly handles GlycoCT and WURCS in this run, but only 27 IUPAC-extended rows converted successfully; GLYCAM-IUPAC, GWB, and compact results from this tool use accession-matched WURCS fallback and are not direct source-parser evidence."
  ),
  collapse = "\n"
)
limitations_body <- paste(
  c(
    "## Limitations and robustness",
    "",
    "Canonical comparison is intentionally limited to chemistry representable by the current glyrepr model. An external output that cannot be normalized is reported as unavailable rather than as a semantic disagreement. Agreement between two direct converters without an accession reference is retained as an unresolved external-consensus difference, not attributed to glyparse, because shared information loss cannot be excluded. The external libraries differ in supported residue vocabularies and reducing-end conventions, so raw string equality is never used.",
    "",
    "Detailed failure classifications reuse exact-row diagnostics only when accession and full source sequence match the prior diagnostic corpus; otherwise the current internal parser stage is replayed. Every final failure row retains the diagnostic source."
  ),
  collapse = "\n"
)
next_steps_body <- paste(
  c(
    "## Recommended next steps",
    "",
    "Review the supported outlier accession(s) first, then the unresolved five-row reference ledger. For parser coverage, prioritize the largest glyrepr-restriction families separately from glyparse implementation limitations. Re-run the workflow after any parser or glyrepr semantic change; the source revisions and corpus hashes in the cache metadata make changes auditable."
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
  path = "data/benchmark/glyparse_validation/format_summary.csv",
  query = list(
    engine = "DuckDB",
    language = "SQL",
    sql = paste0(
      "SELECT * FROM read_csv_auto(",
      "'data/benchmark/glyparse_validation/format_summary.csv')"
    ),
    description = "Read the generated format-level glyparse validation summary."
  )
)
source_tools <- list(
  id = "external_summary",
  label = "External converter summary",
  path = "data/benchmark/glyparse_validation/external_tool_summary.csv",
  query = list(
    engine = "DuckDB",
    language = "SQL",
    sql = paste0(
      "SELECT *, normalized::DOUBLE / rows AS normalized_rate, ",
      "comparable_with_glyparse::DOUBLE / rows AS comparable_rate ",
      "FROM read_csv_auto(",
      "'data/benchmark/glyparse_validation/external_tool_summary.csv')"
    ),
    description = "Read conversion, normalization, and comparison counts by tool and format."
  )
)
source_benchmark <- list(
  id = "parser_benchmark",
  label = "glyparse parser benchmark",
  path = "data/benchmark/glyparse_validation/parser_benchmark.csv",
  query = list(
    engine = "DuckDB",
    language = "SQL",
    sql = paste0(
      "SELECT * FROM read_csv_auto(",
      "'data/benchmark/glyparse_validation/parser_benchmark.csv')"
    ),
    description = "Read median and range statistics from the three parser benchmark repeats."
  )
)
source_failures <- list(
  id = "failure_summary",
  label = "glyparse failure reason summary",
  path = "data/benchmark/glyparse_validation/failure_reason_summary.csv",
  query = list(
    engine = "DuckDB",
    language = "SQL",
    sql = paste0(
      "SELECT * FROM read_csv_auto(",
      "'data/benchmark/glyparse_validation/failure_reason_summary.csv')"
    ),
    description = "Read grouped row-level parser failure classes and reasons."
  )
)
source_disagreements <- list(
  id = "disagreement_ledger",
  label = "Semantic disagreement ledger",
  path = "data/benchmark/glyparse_validation/semantic_disagreements.csv.gz",
  query = list(
    engine = "DuckDB",
    language = "SQL",
    sql = paste0(
      "SELECT * FROM read_csv_auto(",
      "'data/benchmark/glyparse_validation/semantic_disagreements.csv.gz')"
    ),
    description = "Read the complete normalized glyparse-versus-comparator disagreement ledger."
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
          "Canonical equality among accession-matched non-condensed rows."
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
      )
    ),
    sources = list(
      source_summary,
      source_tools,
      source_benchmark,
      source_failures,
      source_disagreements
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
        body = paste(scope_body, methodology_body, sep = "\n\n")
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
        body = paste(limitations_body, next_steps_body, sep = "\n\n")
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
    source_disagreements
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
    "GLYCAM normalization is slower than the other parsers | Single blue root; exact tooltips |"
  )
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
  glycan_format_converter = gfc_metadata
)
jsonlite::write_json(
  run_metadata,
  file.path(benchmark_output_dir, "run-metadata.json"),
  auto_unbox = TRUE,
  pretty = TRUE
)

message("Wrote report artifacts and run-metadata.json")
