#!/usr/bin/env Rscript

script_path <- sub(
  "^--file=",
  "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1L]]
)
source(file.path(dirname(normalizePath(script_path)), "common.R"))

dir.create(benchmark_output_dir, recursive = TRUE, showWarnings = FALSE)

disagreements <- benchmark_read_csv(file.path(
  benchmark_output_dir,
  "semantic_disagreements.csv.gz"
))
external_summary <- benchmark_read_csv(file.path(
  benchmark_output_dir,
  "external_tool_summary.csv"
))
stopifnot(
  nrow(disagreements) > 0L,
  all(nzchar(disagreements$primary_mismatch_class)),
  all(disagreements$mismatch_class_count >= 1L)
)

definitions <- benchmark_mismatch_class_definitions
total_mismatches <- nrow(disagreements)

primary_counts <- table(disagreements$primary_mismatch_class)
reason_summary <- definitions
reason_summary$rows <- as.integer(primary_counts[
  match(reason_summary$mismatch_class, names(primary_counts))
])
reason_summary$rows[is.na(reason_summary$rows)] <- 0L
reason_summary$share <- reason_summary$rows / total_mismatches
reason_summary$unique_accessions <- vapply(
  reason_summary$mismatch_class,
  function(class) {
    length(unique(disagreements$glytoucan_ac[
      disagreements$primary_mismatch_class == class
    ]))
  },
  integer(1)
)
reason_summary$formats <- vapply(
  reason_summary$mismatch_class,
  function(class) {
    length(unique(disagreements$format[
      disagreements$primary_mismatch_class == class
    ]))
  },
  integer(1)
)
reason_summary$comparators <- vapply(
  reason_summary$mismatch_class,
  function(class) {
    length(unique(disagreements$comparator[
      disagreements$primary_mismatch_class == class
    ]))
  },
  integer(1)
)
reason_summary <- reason_summary[order(reason_summary$priority), ]

component_summary <- definitions
component_summary$rows <- vapply(
  component_summary$mismatch_class,
  function(class) {
    sum(disagreements[[paste0("mismatch_has_", class)]])
  },
  integer(1)
)
component_summary$share <- component_summary$rows / total_mismatches
component_summary <- component_summary[
  order(-component_summary$rows, component_summary$priority),
]

tool_reason_summary <- aggregate(
  list(rows = rep(1L, total_mismatches)),
  by = disagreements[, c(
    "comparator",
    "primary_mismatch_class",
    "primary_mismatch_reason"
  )],
  FUN = sum
)
tool_totals <- aggregate(
  list(tool_mismatches = rep(1L, total_mismatches)),
  by = list(comparator = disagreements$comparator),
  FUN = sum
)
tool_reason_summary$tool_mismatches <- tool_totals$tool_mismatches[
  match(tool_reason_summary$comparator, tool_totals$comparator)
]
tool_reason_summary$share_within_tool <-
  tool_reason_summary$rows / tool_reason_summary$tool_mismatches
tool_reason_summary$direct_rows <- vapply(
  seq_len(nrow(tool_reason_summary)),
  function(index) {
    sum(
      disagreements$comparator == tool_reason_summary$comparator[[index]] &
        disagreements$primary_mismatch_class ==
          tool_reason_summary$primary_mismatch_class[[index]] &
        grepl("^direct", disagreements$evidence_scope)
    )
  },
  integer(1)
)
tool_reason_summary$fallback_rows <- vapply(
  seq_len(nrow(tool_reason_summary)),
  function(index) {
    sum(
      disagreements$comparator == tool_reason_summary$comparator[[index]] &
        disagreements$primary_mismatch_class ==
          tool_reason_summary$primary_mismatch_class[[index]] &
        disagreements$evidence_scope == "accession_matched_wurcs_fallback"
    )
  },
  integer(1)
)
tool_reason_summary <- tool_reason_summary[
  order(
    tool_reason_summary$comparator,
    -tool_reason_summary$rows,
    tool_reason_summary$primary_mismatch_class
  ),
]

format_reason_summary <- aggregate(
  list(rows = rep(1L, total_mismatches)),
  by = disagreements[, c(
    "format",
    "primary_mismatch_class",
    "primary_mismatch_reason"
  )],
  FUN = sum
)
format_totals <- aggregate(
  list(format_mismatches = rep(1L, total_mismatches)),
  by = list(format = disagreements$format),
  FUN = sum
)
format_reason_summary$format_mismatches <- format_totals$format_mismatches[
  match(format_reason_summary$format, format_totals$format)
]
format_reason_summary$share_within_format <-
  format_reason_summary$rows / format_reason_summary$format_mismatches
format_reason_summary <- format_reason_summary[
  order(
    format_reason_summary$format,
    -format_reason_summary$rows,
    format_reason_summary$primary_mismatch_class
  ),
]

tool_overall <- aggregate(
  external_summary[, c(
    "comparable_with_glyparse",
    "different_from_glyparse"
  )],
  by = list(comparator = external_summary$tool),
  FUN = sum
)
tool_overall$mismatch_rate <-
  tool_overall$different_from_glyparse /
  tool_overall$comparable_with_glyparse
tool_overall <- tool_overall[order(-tool_overall$different_from_glyparse), ]

example_order <- order(
  match(disagreements$primary_mismatch_class, definitions$mismatch_class),
  disagreements$comparator,
  disagreements$format,
  disagreements$row_index
)
ordered_disagreements <- disagreements[example_order, ]
example_rows <- ordered_disagreements[
  !duplicated(ordered_disagreements$primary_mismatch_class),
  c(
    "primary_mismatch_class",
    "primary_mismatch_reason",
    "format",
    "comparator",
    "glytoucan_ac",
    "difference_reason"
  )
]
example_rows <- example_rows[
  order(match(example_rows$primary_mismatch_class, definitions$mismatch_class)),
]

multi_component_rows <- sum(disagreements$mismatch_class_count > 1L)
largest <- reason_summary[which.max(reason_summary$rows), ]
unique_accessions <- length(unique(disagreements$glytoucan_ac))
headline <- data.frame(
  total_mismatches = total_mismatches,
  primary_classes = nrow(reason_summary),
  largest_primary_share = largest$share,
  multi_component_share = multi_component_rows / total_mismatches,
  unique_accessions = unique_accessions,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

benchmark_write_csv(
  reason_summary,
  file.path(
    benchmark_output_dir,
    "semantic_disagreement_reason_summary.csv"
  )
)
benchmark_write_csv(
  component_summary,
  file.path(
    benchmark_output_dir,
    "semantic_disagreement_component_summary.csv"
  )
)
benchmark_write_csv(
  tool_reason_summary,
  file.path(
    benchmark_output_dir,
    "semantic_disagreement_tool_reason_summary.csv"
  )
)
benchmark_write_csv(
  format_reason_summary,
  file.path(
    benchmark_output_dir,
    "semantic_disagreement_format_reason_summary.csv"
  )
)
benchmark_write_csv(
  example_rows,
  file.path(
    benchmark_output_dir,
    "semantic_disagreement_examples.csv"
  )
)

format_integer <- function(value) {
  format(value, big.mark = ",", scientific = FALSE)
}

tool_top_rows <- do.call(
  rbind,
  lapply(
    split(tool_reason_summary, tool_reason_summary$comparator),
    function(x) {
      x[which.max(x$rows), ]
    }
  )
)
rownames(tool_top_rows) <- NULL
tool_lines <- vapply(
  seq_len(nrow(tool_top_rows)),
  function(index) {
    row <- tool_top_rows[index, ]
    paste0(
      "- **",
      row$comparator,
      ":** ",
      row$primary_mismatch_reason,
      " is the largest class at ",
      format_integer(row$rows),
      " rows (",
      sprintf("%.1f%%", 100 * row$share_within_tool),
      " of this tool's mismatches)."
    )
  },
  character(1)
)

generated_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
technical_summary_body <- paste(
  c(
    "## Technical summary",
    "",
    paste0(
      "The classified ledger contains **",
      format_integer(total_mismatches),
      "** normalized glyparse-versus-tool mismatches across **",
      format_integer(unique_accessions),
      "** GlyTouCan accessions. Every row has one primary reason and one or ",
      "more component flags; no row is unclassified."
    ),
    "",
    paste0(
      "The largest primary class is **",
      largest$mismatch_reason,
      "** with **",
      format_integer(largest$rows),
      "** rows (**",
      sprintf("%.1f%%", 100 * largest$share),
      "**). **",
      format_integer(multi_component_rows),
      "** rows (**",
      sprintf("%.1f%%", 100 * multi_component_rows / total_mismatches),
      "**) carry more than one component flag, so the overlapping component ",
      "totals intentionally exceed the ledger row count."
    )
  ),
  collapse = "\n"
)
definition_body <- paste(
  c(
    "## Classification contract",
    "",
    "Primary reasons are mutually exclusive and follow a specificity-first priority: reducing-end alditol, floating structure/domain, reducing-end anomer, residue composition, topology/attachment, internal linkage/anomer, then same-shape label assignment. This ordering makes the primary totals reconcile exactly while preserving every detected difference in `mismatch_classes` and seven boolean component fields.",
    "",
    "The classes describe how the two normalized glycan graphs differ. They do not by themselves identify which tool is correct or whether the difference reflects expected information loss, a source-format ambiguity, or a converter defect."
  ),
  collapse = "\n"
)
primary_body <- paste(
  c(
    "## Reducing-end and composition differences dominate",
    "",
    paste0(
      "The ranked chart partitions all ",
      format_integer(total_mismatches),
      " rows into exactly one primary class. Reducing-end alditol, internal ",
      "linkage/anomer, and residue/substituent composition together account ",
      "for **",
      sprintf(
        "%.1f%%",
        100 *
          sum(reason_summary$share[
            reason_summary$mismatch_class %in%
              c(
                "reducing_end_alditol",
                "linkage_or_anomer",
                "residue_composition"
              )
          ])
      ),
      "** of mismatches. This concentrates manual review on information-loss ",
      "and chemistry-normalization behavior before rarer topology cases."
    )
  ),
  collapse = "\n"
)
tool_body <- paste(
  c(
    "## The mismatch mix is tool-specific",
    "",
    "The grouped chart uses the same primary partition for each comparator; its bars should be read as counts, not conversion coverage. Direct and accession-matched fallback counts remain separate in the exact table.",
    "",
    tool_lines
  ),
  collapse = "\n"
)
component_body <- paste(
  c(
    "## Component flags preserve compound differences",
    "",
    paste0(
      "A row may differ in composition, linkage, alditol state, topology, and ",
      "floating metadata simultaneously. The component table therefore uses ",
      "overlapping incidence counts: **",
      format_integer(component_summary$rows[[1L]]),
      "** rows include the most common component, **",
      component_summary$mismatch_reason[[1L]],
      "**. These counts diagnose what changed without forcing a single-cause ",
      "interpretation."
    )
  ),
  collapse = "\n"
)
examples_body <- paste(
  c(
    "## Representative rows make each class auditable",
    "",
    "The example table retains the first deterministic row for each primary class after sorting by class priority, comparator, format, and source row. It is an index into the complete compressed ledger, not a substitute for row-level review."
  ),
  collapse = "\n"
)
limitations_body <- paste(
  c(
    "## Limitations and robustness",
    "",
    "Classification is deterministic over the component-level reasons generated by the graph comparator, and the workflow fails if any row matches no known class. Mismatches are included only when both glyparse and the external output normalized successfully; conversion failures and outputs outside the current glyrepr model are excluded. Counts include both direct source conversions and explicitly labeled accession-matched WURCS fallbacks, so they should not be interpreted as direct-parser defect rates without filtering `evidence_scope`."
  ),
  collapse = "\n"
)
next_steps_body <- paste(
  c(
    "## Recommended next steps",
    "",
    "- Review reducing-end alditol differences in glycowork separately from ordinary residue-composition differences.",
    "- Review glypy residue/substituent composition differences for vocabulary or serialization loss.",
    "- Review GlycanFormatConverter reducing-end-anomer differences with the strict contract and keep fallback evidence separate.",
    "- Manually inspect the rare topology/attachment and same-shape label-assignment classes first because they are more likely to reveal structural conversion defects."
  ),
  collapse = "\n"
)
questions_body <- paste(
  c(
    "## Further questions",
    "",
    "Which primary classes should be treated as expected tool information loss versus defects requiring upstream reports? That attribution needs format-aware manual adjudication or an additional independent reference; the current classification intentionally stops at graph-difference type."
  ),
  collapse = "\n"
)

report_markdown <- c(
  "# Semantic disagreement classification",
  "",
  technical_summary_body,
  "",
  definition_body,
  "",
  primary_body,
  "",
  tool_body,
  "",
  component_body,
  "",
  examples_body,
  "",
  limitations_body,
  "",
  next_steps_body,
  "",
  questions_body
)
writeLines(
  report_markdown,
  file.path(
    benchmark_output_dir,
    "semantic_disagreements_report.md"
  ),
  useBytes = TRUE
)

source_ledger <- list(
  id = "disagreement_ledger",
  label = "Classified semantic disagreement ledger",
  path = "results/data/glyparse_validation/semantic_disagreements.csv.gz",
  query = list(
    engine = "DuckDB",
    language = "SQL",
    sql = paste0(
      "SELECT * FROM read_csv_auto(",
      "'results/data/glyparse_validation/semantic_disagreements.csv.gz')"
    ),
    description = "Read all classified normalized glyparse-versus-tool mismatches."
  )
)
source_reason <- list(
  id = "reason_summary",
  label = "Primary semantic mismatch reason summary",
  path = paste0(
    "results/data/glyparse_validation/",
    "semantic_disagreement_reason_summary.csv"
  ),
  query = list(
    engine = "DuckDB",
    language = "SQL",
    sql = paste0(
      "SELECT * FROM read_csv_auto(",
      "'results/data/glyparse_validation/",
      "semantic_disagreement_reason_summary.csv')"
    ),
    description = "Read mutually exclusive primary reason counts and shares."
  )
)
source_tool <- list(
  id = "tool_reason_summary",
  label = "Tool by primary mismatch reason summary",
  path = paste0(
    "results/data/glyparse_validation/",
    "semantic_disagreement_tool_reason_summary.csv"
  ),
  query = list(
    engine = "DuckDB",
    language = "SQL",
    sql = paste0(
      "SELECT * FROM read_csv_auto(",
      "'results/data/glyparse_validation/",
      "semantic_disagreement_tool_reason_summary.csv')"
    ),
    description = "Read primary mismatch counts by comparator and evidence scope."
  )
)
source_component <- list(
  id = "component_summary",
  label = "Overlapping semantic mismatch component summary",
  path = paste0(
    "results/data/glyparse_validation/",
    "semantic_disagreement_component_summary.csv"
  ),
  query = list(
    engine = "DuckDB",
    language = "SQL",
    sql = paste0(
      "SELECT * FROM read_csv_auto(",
      "'results/data/glyparse_validation/",
      "semantic_disagreement_component_summary.csv')"
    ),
    description = "Read overlapping mismatch component incidence counts."
  )
)
source_external <- list(
  id = "external_summary",
  label = "External comparison summary",
  path = "results/data/glyparse_validation/external_tool_summary.csv",
  query = list(
    engine = "DuckDB",
    language = "SQL",
    sql = paste0(
      "SELECT * FROM read_csv_auto(",
      "'results/data/glyparse_validation/external_tool_summary.csv')"
    ),
    description = "Read comparable and different normalized rows by tool and format."
  )
)
sources <- list(
  source_ledger,
  source_reason,
  source_tool,
  source_component,
  source_external
)

artifact <- list(
  surface = "report",
  manifest = list(
    version = 1L,
    surface = "report",
    title = "Semantic disagreement classification",
    description = paste0(
      "Reason classification for normalized glyparse-versus-tool semantic ",
      "mismatches."
    ),
    generatedAt = generated_at,
    cards = list(
      list(
        id = "mismatches",
        description = "All normalized glyparse-versus-tool differences.",
        dataset = "headline",
        sourceId = "disagreement_ledger",
        metrics = list(list(
          label = "Classified mismatches",
          field = "total_mismatches",
          format = "number"
        ))
      ),
      list(
        id = "classes",
        description = "Mutually exclusive primary mismatch classes.",
        dataset = "headline",
        sourceId = "reason_summary",
        metrics = list(list(
          label = "Primary classes",
          field = "primary_classes",
          format = "number"
        ))
      ),
      list(
        id = "largest",
        description = paste0(
          "Share assigned to the largest primary class: ",
          largest$mismatch_reason,
          "."
        ),
        dataset = "headline",
        sourceId = "reason_summary",
        metrics = list(list(
          label = "Largest class share",
          field = "largest_primary_share",
          format = "percent"
        ))
      ),
      list(
        id = "compound",
        description = "Rows carrying more than one component-level mismatch flag.",
        dataset = "headline",
        sourceId = "disagreement_ledger",
        metrics = list(list(
          label = "Multi-component share",
          field = "multi_component_share",
          format = "percent"
        ))
      )
    ),
    charts = list(
      list(
        id = "primary_reasons",
        title = "Primary mismatch reasons",
        subtitle = paste0(
          "Mutually exclusive partition of ",
          format_integer(total_mismatches),
          " normalized mismatches."
        ),
        type = "bar",
        dataset = "reason_summary",
        sourceId = "reason_summary",
        encodings = list(
          x = list(
            field = "mismatch_reason",
            type = "nominal",
            label = "Primary reason"
          ),
          y = list(
            field = "rows",
            type = "quantitative",
            label = "Mismatch rows",
            format = "number"
          )
        )
      ),
      list(
        id = "tool_reasons",
        title = "Primary mismatch reasons by comparator",
        subtitle = paste0(
          "Counts among normalized comparable rows; direct and fallback scope ",
          "is retained in the supporting table."
        ),
        type = "bar",
        dataset = "tool_reason_summary",
        sourceId = "tool_reason_summary",
        encodings = list(
          x = list(
            field = "primary_mismatch_reason",
            type = "nominal",
            label = "Primary reason"
          ),
          y = list(
            field = "rows",
            type = "quantitative",
            label = "Mismatch rows",
            format = "number"
          ),
          color = list(
            field = "comparator",
            type = "nominal",
            label = "Comparator"
          )
        )
      )
    ),
    tables = list(
      list(
        id = "reason_detail",
        title = "Primary reason definitions and counts",
        dataset = "reason_summary",
        sourceId = "reason_summary",
        defaultSort = list(field = "rows", direction = "desc"),
        columns = list(
          list(
            field = "mismatch_reason",
            label = "Primary reason",
            type = "text"
          ),
          list(field = "definition", label = "Definition", type = "text"),
          list(field = "rows", label = "Rows", format = "number"),
          list(field = "share", label = "Share", format = "percent"),
          list(
            field = "unique_accessions",
            label = "Accessions",
            format = "number"
          ),
          list(field = "formats", label = "Formats", format = "number"),
          list(field = "comparators", label = "Tools", format = "number")
        )
      ),
      list(
        id = "tool_detail",
        title = "Tool-specific primary reason mix",
        dataset = "tool_reason_summary",
        sourceId = "tool_reason_summary",
        defaultSort = list(field = "rows", direction = "desc"),
        columns = list(
          list(field = "comparator", label = "Comparator", type = "text"),
          list(
            field = "primary_mismatch_reason",
            label = "Primary reason",
            type = "text"
          ),
          list(field = "rows", label = "Rows", format = "number"),
          list(
            field = "share_within_tool",
            label = "Within-tool share",
            format = "percent"
          ),
          list(field = "direct_rows", label = "Direct", format = "number"),
          list(field = "fallback_rows", label = "Fallback", format = "number")
        )
      ),
      list(
        id = "component_detail",
        title = "Overlapping component flags",
        dataset = "component_summary",
        sourceId = "component_summary",
        defaultSort = list(field = "rows", direction = "desc"),
        columns = list(
          list(field = "mismatch_reason", label = "Component", type = "text"),
          list(field = "definition", label = "Definition", type = "text"),
          list(field = "rows", label = "Rows", format = "number"),
          list(field = "share", label = "Incidence", format = "percent")
        )
      ),
      list(
        id = "examples",
        title = "Representative mismatch rows",
        dataset = "examples",
        sourceId = "disagreement_ledger",
        defaultSort = list(
          field = "primary_mismatch_reason",
          direction = "asc"
        ),
        columns = list(
          list(
            field = "primary_mismatch_reason",
            label = "Primary reason",
            type = "text"
          ),
          list(field = "format", label = "Format", type = "text"),
          list(field = "comparator", label = "Comparator", type = "text"),
          list(field = "glytoucan_ac", label = "GlyTouCan", type = "text"),
          list(
            field = "difference_reason",
            label = "Detailed component reason",
            type = "text"
          )
        )
      )
    ),
    sources = sources,
    blocks = list(
      list(
        id = "title",
        type = "markdown",
        body = "# Semantic disagreement classification"
      ),
      list(
        id = "technical_summary",
        type = "markdown",
        body = technical_summary_body,
        sourceId = "disagreement_ledger"
      ),
      list(
        id = "headline_metrics",
        type = "metric-strip",
        cardIds = c("mismatches", "classes", "largest", "compound")
      ),
      list(
        id = "classification_contract",
        type = "markdown",
        body = definition_body
      ),
      list(
        id = "primary_findings",
        type = "markdown",
        body = primary_body,
        sourceId = "reason_summary"
      ),
      list(
        id = "primary_chart",
        type = "chart",
        chartId = "primary_reasons",
        layout = "full"
      ),
      list(
        id = "reason_table",
        type = "table",
        tableId = "reason_detail"
      ),
      list(
        id = "tool_findings",
        type = "markdown",
        body = tool_body,
        sourceId = "tool_reason_summary"
      ),
      list(
        id = "tool_chart",
        type = "chart",
        chartId = "tool_reasons",
        layout = "full"
      ),
      list(
        id = "tool_table",
        type = "table",
        tableId = "tool_detail"
      ),
      list(
        id = "component_findings",
        type = "markdown",
        body = component_body,
        sourceId = "component_summary"
      ),
      list(
        id = "component_table",
        type = "table",
        tableId = "component_detail"
      ),
      list(
        id = "examples_findings",
        type = "markdown",
        body = examples_body
      ),
      list(
        id = "examples_table",
        type = "table",
        tableId = "examples"
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
      ),
      list(
        id = "further_questions",
        type = "markdown",
        body = questions_body
      )
    )
  ),
  snapshot = list(
    version = 1L,
    generatedAt = generated_at,
    status = "ready",
    datasets = list(
      headline = headline,
      reason_summary = reason_summary,
      tool_reason_summary = tool_reason_summary,
      component_summary = component_summary,
      tool_overall = tool_overall,
      examples = example_rows
    )
  ),
  sources = sources
)

jsonlite::write_json(
  artifact,
  file.path(
    benchmark_output_dir,
    "semantic_disagreements_artifact.json"
  ),
  auto_unbox = TRUE,
  pretty = TRUE,
  na = "null",
  digits = 12
)

chart_map <- c(
  "# Semantic disagreement chart map",
  "",
  paste0(
    "| Segment | Question | Family/type | Fields | Supported claim | ",
    "Palette |"
  ),
  "|---|---|---|---|---|---|",
  paste0(
    "| Primary reasons | How are all mismatch rows partitioned? | ",
    "Comparison / ranked bar | mismatch_reason, rows, share | ",
    "Reducing-end and composition reasons dominate | Single blue root; ",
    "direct axis labels |"
  ),
  paste0(
    "| Tool mix | How does the primary reason mix differ by comparator? | ",
    "Comparison / grouped bar | primary_mismatch_reason, rows, ",
    "comparator, direct_rows, fallback_rows | Each tool has a distinct ",
    "dominant class | Three approved roots plus legend |"
  ),
  "",
  paste0(
    "Overlapping component incidence, exact reason definitions, evidence ",
    "scope, and representative rows are tables because auditability matters ",
    "more than another repeated bar chart."
  )
)
writeLines(
  chart_map,
  file.path(
    benchmark_output_dir,
    "semantic_disagreements_chart_map.md"
  ),
  useBytes = TRUE
)

stopifnot(
  sum(reason_summary$rows) == total_mismatches,
  all(reason_summary$rows > 0L),
  all(
    component_summary$rows >=
      reason_summary$rows[
        match(
          component_summary$mismatch_class,
          reason_summary$mismatch_class
        )
      ]
  ),
  sum(tool_reason_summary$rows) == total_mismatches,
  sum(format_reason_summary$rows) == total_mismatches,
  nrow(example_rows) == nrow(reason_summary),
  all(nzchar(example_rows$difference_reason))
)

message("Wrote semantic disagreement classification report artifacts.")
