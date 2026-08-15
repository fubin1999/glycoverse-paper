benchmark_script_dir <- function() {
  command <- commandArgs(trailingOnly = FALSE)
  file_argument <- grep("^--file=", command, value = TRUE)
  if (length(file_argument) == 0L) {
    return(normalizePath(getwd(), mustWork = TRUE))
  }
  dirname(normalizePath(sub("^--file=", "", file_argument[[1L]])))
}

benchmark_module_dir <- benchmark_script_dir()
benchmark_project_dir <- normalizePath(
  Sys.getenv(
    "GLYPARSE_BENCHMARK_ROOT",
    unset = file.path(benchmark_module_dir, "..", "..", "..")
  ),
  mustWork = TRUE
)
benchmark_corpus_dir <- file.path(
  benchmark_project_dir,
  "data",
  "benchmark",
  "glycan_sequences"
)
benchmark_output_dir <- file.path(
  benchmark_project_dir,
  "results",
  "data",
  "glyparse_validation"
)
benchmark_cache_dir <- file.path(benchmark_output_dir, "cache")
benchmark_cores <- suppressWarnings(as.integer(
  Sys.getenv("GLYPARSE_BENCHMARK_CORES", unset = "8")
))
if (is.na(benchmark_cores) || benchmark_cores < 1L) {
  stop("GLYPARSE_BENCHMARK_CORES must be a positive integer")
}

benchmark_glyparse_dir <- normalizePath(
  Sys.getenv(
    "GLYPARSE_SOURCE_DIR",
    unset = "/Users/fubin/glycoverse/packages/glyparse"
  ),
  mustWork = TRUE
)
benchmark_glyrepr_dir <- normalizePath(
  Sys.getenv(
    "GLYREPR_SOURCE_DIR",
    unset = "/Users/fubin/glycoverse/packages/glyrepr"
  ),
  mustWork = TRUE
)

benchmark_definitions <- list(
  glycam_iupac = list(
    file = "glycan_sequences_glycam_iupac.csv",
    sequence_column = "sequence_glycam_iupac",
    parser = "parse_glycam_iupac",
    label = "GLYCAM-IUPAC",
    reference_validation = TRUE
  ),
  glycoct = list(
    file = "glycan_sequences_glycoct.csv",
    sequence_column = "sequence_glycoct",
    parser = "parse_glycoct",
    label = "GlycoCT",
    reference_validation = TRUE
  ),
  gwb = list(
    file = "glycan_sequences_gwb.csv",
    sequence_column = "sequence_gwb",
    parser = "parse_gwb",
    label = "GlycoWorkbench",
    reference_validation = TRUE
  ),
  iupac_compact = list(
    file = "glycan_sequences_iupac_compact.csv",
    sequence_column = "sequence_iupac_compact",
    parser = "parse_iupac_compact",
    label = "IUPAC-compact",
    reference_validation = TRUE
  ),
  iupac_condensed = list(
    file = "glycan_sequences_iupac_condensed.csv",
    sequence_column = "sequence_iupac_condensed",
    parser = "parse_iupac_condensed",
    label = "IUPAC-condensed",
    reference_validation = FALSE
  ),
  iupac_extended = list(
    file = "glycan_sequences_iupac_extended.csv",
    sequence_column = "sequence_iupac_extended",
    parser = "parse_iupac_extended",
    label = "IUPAC-extended",
    reference_validation = TRUE
  ),
  wurcs = list(
    file = "glycan_sequences_wurcs.csv",
    sequence_column = "sequence_wurcs",
    parser = "parse_wurcs",
    label = "WURCS",
    reference_validation = TRUE
  )
)

benchmark_load_packages <- function() {
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("The pkgload package is required")
  }
  pkgload::load_all(benchmark_glyrepr_dir, reset = TRUE, quiet = TRUE)
  pkgload::load_all(benchmark_glyparse_dir, reset = TRUE, quiet = TRUE)
  stopifnot(
    identical(
      normalizePath(getNamespaceInfo(asNamespace("glyrepr"), "path")),
      benchmark_glyrepr_dir
    ),
    identical(
      normalizePath(getNamespaceInfo(asNamespace("glyparse"), "path")),
      benchmark_glyparse_dir
    )
  )
}

benchmark_revision <- function(path) {
  output <- system2(
    "git",
    c("-C", path, "rev-parse", "HEAD"),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) {
    return("unavailable")
  }
  trimws(output[[1L]])
}

benchmark_sha256 <- function(path) {
  output <- system2(
    "shasum",
    c("-a", "256", path),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) {
    stop("Unable to hash ", path)
  }
  strsplit(output[[1L]], "[[:space:]]+")[[1L]][[1L]]
}

benchmark_read_corpus <- function(definition) {
  path <- file.path(benchmark_corpus_dir, definition$file)
  data <- utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    comment.char = ""
  )
  expected <- c("glytoucan_ac", definition$sequence_column)
  stopifnot(
    identical(names(data), expected),
    !anyNA(data),
    !anyDuplicated(data$glytoucan_ac),
    all(nzchar(data$glytoucan_ac)),
    all(nzchar(data[[definition$sequence_column]]))
  )
  data
}

benchmark_preprocess <- function(format, sequence) {
  if (identical(format, "glycoct")) {
    return(gsub(" ", "\n", sequence, fixed = TRUE))
  }
  sequence
}

benchmark_condition_text <- function(condition) {
  text <- conditionMessage(condition)
  text <- gsub("\033\\[[0-9;]*m", "", text)
  text <- gsub("[\r\n]+", " ", text)
  trimws(gsub("[[:space:]]+", " ", text))
}

benchmark_parser_call <- function(parser, sequence) {
  suppressWarnings(do.call(
    parser,
    list(x = sequence, on_failure = "na", progress = FALSE)
  ))
}

benchmark_parse_recursive <- function(parser, sequence) {
  output <- tryCatch(
    benchmark_parser_call(parser, sequence),
    error = identity
  )
  if (!inherits(output, "condition")) {
    canonical <- as.character(output)
    return(data.frame(
      parse_status = ifelse(is.na(canonical), "failed", "parsed"),
      parsed_iupac_condensed = ifelse(is.na(canonical), "", canonical),
      parse_error = "",
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }
  if (length(sequence) == 1L) {
    return(data.frame(
      parse_status = "failed",
      parsed_iupac_condensed = "",
      parse_error = benchmark_condition_text(output),
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }
  middle <- length(sequence) %/% 2L
  rbind(
    benchmark_parse_recursive(parser, sequence[seq_len(middle)]),
    benchmark_parse_recursive(parser, sequence[(middle + 1L):length(sequence)])
  )
}

benchmark_capture_failure <- function(parser, sequence) {
  tryCatch(
    {
      suppressWarnings(do.call(
        parser,
        list(x = sequence, on_failure = "error", progress = FALSE)
      ))
      "Parser returned a missing value without exposing an underlying condition."
    },
    error = benchmark_condition_text
  )
}

benchmark_parse_one <- function(parser, sequence) {
  tryCatch(
    {
      output <- suppressWarnings(do.call(
        parser,
        list(x = sequence, on_failure = "error", progress = FALSE)
      ))
      canonical <- as.character(output)
      if (length(canonical) != 1L || is.na(canonical)) {
        stop("Parser returned a missing or non-scalar result")
      }
      data.frame(
        parse_status = "parsed",
        parsed_iupac_condensed = canonical,
        mono_type = glyrepr::get_mono_type(output),
        parse_error = "",
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    },
    error = function(condition) {
      data.frame(
        parse_status = "failed",
        parsed_iupac_condensed = "",
        mono_type = "",
        parse_error = benchmark_condition_text(condition),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  )
}

benchmark_parse_corpus <- function(
  format,
  definition,
  data,
  chunk_size = 1000L
) {
  parser <- getExportedValue("glyparse", definition$parser)
  sequence <- benchmark_preprocess(
    format,
    data[[definition$sequence_column]]
  )
  chunks <- split(
    seq_along(sequence),
    ceiling(seq_along(sequence) / chunk_size)
  )
  results <- vector("list", length(chunks))
  started <- proc.time()[["elapsed"]]
  for (index in seq_along(chunks)) {
    rows <- chunks[[index]]
    scalar_results <- parallel::mclapply(
      sequence[rows],
      function(value) benchmark_parse_one(parser, value),
      mc.cores = benchmark_cores,
      mc.preschedule = TRUE
    )
    results[[index]] <- do.call(rbind, scalar_results)
    message(
      definition$label,
      ": ",
      max(rows),
      "/",
      length(sequence),
      " rows parsed"
    )
  }
  result <- do.call(rbind, results)
  elapsed <- proc.time()[["elapsed"]] - started
  list(result = result, elapsed = elapsed)
}

benchmark_sample_indices <- function(size, target = 1000L) {
  if (size <= target) {
    return(seq_len(size))
  }
  unique(as.integer(round(seq(1, size, length.out = target))))
}

benchmark_measure_parser <- function(
  format,
  definition,
  data,
  validation_result,
  repeats = 3L
) {
  parser <- getExportedValue("glyparse", definition$parser)
  parsed_rows <- validation_result$parse_status == "parsed"
  type_counts <- sort(
    table(validation_result$mono_type[parsed_rows]),
    decreasing = TRUE
  )
  benchmark_type <- names(type_counts)[[1L]]
  candidates <- which(
    parsed_rows & validation_result$mono_type == benchmark_type
  )
  indices <- candidates[benchmark_sample_indices(length(candidates))]
  sequence <- benchmark_preprocess(
    format,
    data[[definition$sequence_column]][indices]
  )
  elapsed <- numeric(repeats)
  parsed <- integer(repeats)
  for (index in seq_len(repeats)) {
    gc()
    started <- proc.time()[["elapsed"]]
    result <- benchmark_parser_call(parser, sequence)
    elapsed[[index]] <- proc.time()[["elapsed"]] - started
    parsed[[index]] <- sum(!is.na(result))
  }
  data.frame(
    format = format,
    parser = definition$parser,
    mono_type = benchmark_type,
    sample_rows = length(sequence),
    repeat_id = seq_len(repeats),
    elapsed_seconds = elapsed,
    rows_per_second = length(sequence) / elapsed,
    parsed_rows = parsed,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

benchmark_write_csv <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  connection <- if (grepl("\\.gz$", path)) {
    gzfile(path, open = "wt", encoding = "UTF-8")
  } else {
    file(path, open = "wt", encoding = "UTF-8")
  }
  on.exit(close(connection), add = TRUE)
  utils::write.csv(
    data,
    connection,
    row.names = FALSE,
    na = "",
    quote = TRUE
  )
}

benchmark_read_csv <- function(path) {
  utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    comment.char = ""
  )
}

benchmark_parse_iupac_values <- function(values) {
  parser <- getExportedValue("glyparse", "parse_iupac_condensed")
  scalar_results <- parallel::mclapply(
    values,
    function(value) benchmark_parse_one(parser, value),
    mc.cores = benchmark_cores,
    mc.preschedule = TRUE
  )
  do.call(rbind, scalar_results)
}

benchmark_normalize_iupac <- function(values, chunk_size = 1000L) {
  result <- data.frame(
    normalization_status = rep("not_attempted", length(values)),
    normalized_iupac_condensed = rep("", length(values)),
    normalization_error = rep("", length(values)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  available <- which(!is.na(values) & nzchar(values))
  if (length(available) == 0L) {
    return(result)
  }
  chunks <- split(
    available,
    ceiling(seq_along(available) / chunk_size)
  )
  for (rows in chunks) {
    parsed <- benchmark_parse_iupac_values(values[rows])
    ok <- parsed$parse_status == "parsed"
    result$normalization_status[rows] <- ifelse(ok, "normalized", "failed")
    result$normalized_iupac_condensed[rows] <- parsed$parsed_iupac_condensed
    result$normalization_error[rows] <- parsed$parse_error
  }
  result
}

benchmark_normalize_iupac_fast <- function(values, chunk_size = 500L) {
  result <- data.frame(
    normalization_status = rep("not_attempted", length(values)),
    normalized_iupac_condensed = rep("", length(values)),
    normalization_error = rep("", length(values)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  available <- which(!is.na(values) & nzchar(values))
  if (length(available) == 0L) {
    return(result)
  }
  generic_pattern <- paste0(
    "(?<![[:alnum:]])",
    "(?:HexNAc|HexN|HexA|Hex|dHex|Pen|Sia|Unk)",
    "(?=[0-9\\(\\{\\[?]|$)"
  )
  generic <- grepl(generic_pattern, values[available], perl = TRUE)
  groups <- list(available[!generic], available[generic])
  chunks <- unlist(
    lapply(groups, function(rows) {
      if (length(rows) == 0L) {
        return(list())
      }
      split(rows, ceiling(seq_along(rows) / chunk_size))
    }),
    recursive = FALSE
  )
  parser <- getExportedValue("glyparse", "parse_iupac_condensed")
  parsed_chunks <- parallel::mclapply(
    chunks,
    function(rows) {
      list(
        rows = rows,
        parsed = benchmark_parse_recursive(parser, values[rows])
      )
    },
    mc.cores = benchmark_cores,
    mc.preschedule = TRUE
  )
  completed <- 0L
  for (chunk in parsed_chunks) {
    rows <- chunk$rows
    parsed <- chunk$parsed
    ok <- parsed$parse_status == "parsed"
    result$normalization_status[rows] <- ifelse(ok, "normalized", "failed")
    result$normalized_iupac_condensed[rows] <- parsed$parsed_iupac_condensed
    result$normalization_error[rows] <- parsed$parse_error
    completed <- completed + length(rows)
    if (completed %% 10000L < length(rows)) {
      message("External normalization: ", completed, "/", length(available))
    }
  }
  failed <- which(
    result$normalization_status == "failed" &
      !nzchar(result$normalization_error)
  )
  if (length(failed) > 0L) {
    result$normalization_error[failed] <- paste0(
      "parse_iupac_condensed() returned NA during vector normalization; ",
      "the converter output is not comparable in the current canonical model."
    )
  }
  result
}

benchmark_strip_reducing_end_token <- function(value) {
  sub("\\([ab?][0-9?/]+-$", "", value, perl = TRUE)
}

benchmark_comparison_key <- function(value, comparison_contract) {
  if (
    identical(
      comparison_contract,
      "graph_isomorphism_ignore_terminal_reducing_end_token"
    )
  ) {
    return(benchmark_strip_reducing_end_token(value))
  }
  value
}

benchmark_semantic_inventory_key <- function(value, comparison_contract) {
  residues <- regmatches(
    value,
    gregexpr("[[:alnum:]?+-]+(?=\\()", value, perl = TRUE)
  )[[1L]]
  linkages <- regmatches(
    value,
    gregexpr("(?<=\\()[ab?][0-9?/]+-[0-9?/]*", value, perl = TRUE)
  )[[1L]]
  if (identical(residues, "")) {
    residues <- character()
  }
  if (identical(linkages, "")) {
    linkages <- character()
  }
  if (
    identical(
      comparison_contract,
      "graph_isomorphism_ignore_terminal_reducing_end_token"
    ) &&
      grepl("\\([ab?][0-9?/]+-$", value, perl = TRUE) &&
      length(linkages) > 0L
  ) {
    linkages <- head(linkages, -1L)
  }
  linkages <- sub(
    "^([ab?])[0-9?/]+-",
    "\\1*-",
    linkages,
    perl = TRUE
  )
  paste(
    paste(sort(residues, method = "radix"), collapse = "\r"),
    paste(sort(linkages, method = "radix"), collapse = "\r"),
    grepl("-ol\\(", value),
    sep = "\n"
  )
}

benchmark_semantic_fingerprint_cache <- new.env(parent = emptyenv())

benchmark_semantic_fingerprint <- function(value, comparison_contract) {
  cache_key <- paste(comparison_contract, value, sep = "\r")
  if (
    exists(
      cache_key,
      envir = benchmark_semantic_fingerprint_cache,
      inherits = FALSE
    )
  ) {
    return(get(
      cache_key,
      envir = benchmark_semantic_fingerprint_cache,
      inherits = FALSE
    ))
  }

  fingerprint <- tryCatch(
    {
      structure <- glyparse::parse_iupac_condensed(
        value,
        on_failure = "error",
        progress = FALSE
      )
      structure <- glyrepr::fill_anomer_pos(structure)
      graph <- glyrepr::get_structure_graphs(structure)
      parts <- glyrepr:::normalize_floating_parts(graph)
      substituents <- glyrepr:::normalize_floating_substituents(graph)
      main_vertices <- glyrepr:::floating_metadata_main_vertices(graph, parts)
      labels <- glyrepr:::floating_augmented_structure_labels(
        graph,
        parts,
        substituents,
        main_vertices
      )
      vertex_labels <- labels$vertices
      endpoints <- igraph::as_edgelist(graph, names = FALSE)
      anomer <- if (
        identical(
          comparison_contract,
          "graph_isomorphism_ignore_terminal_reducing_end_token"
        )
      ) {
        "ignored"
      } else {
        igraph::graph_attr(graph, "anomer")
      }
      records <- c(
        paste(
          "graph",
          anomer,
          isTRUE(igraph::graph_attr(graph, "alditol")),
          sep = "\r"
        ),
        paste(
          "residue",
          vertex_labels,
          igraph::V(graph)$mono,
          igraph::V(graph)$sub,
          sep = "\r"
        )
      )
      if (nrow(endpoints) > 0L) {
        records <- c(
          records,
          paste(
            "edge",
            vertex_labels[endpoints[, 1L]],
            vertex_labels[endpoints[, 2L]],
            igraph::E(graph)$linkage,
            sep = "\r"
          )
        )
      }
      for (index in seq_along(parts)) {
        part <- parts[[index]]
        candidates <- glyrepr:::floating_part_candidate_parents(graph, part)
        records <- c(
          records,
          paste(
            "floating-part",
            labels$parts[[index]],
            vertex_labels[[part$root]],
            part$linkage,
            paste(
              sort(vertex_labels[candidates], method = "radix"),
              collapse = ","
            ),
            sep = "\r"
          )
        )
      }
      for (index in seq_along(substituents)) {
        substituent <- substituents[[index]]
        candidates <- glyrepr:::floating_substituent_candidate_parents(
          graph,
          substituent
        )
        records <- c(
          records,
          paste(
            "floating-substituent",
            labels$substituents[[index]],
            substituent$substituent,
            paste(
              sort(vertex_labels[candidates], method = "radix"),
              collapse = ","
            ),
            sep = "\r"
          )
        )
      }
      paste(sort(records, method = "radix"), collapse = "\n")
    },
    error = function(condition) NA_character_
  )
  assign(
    cache_key,
    fingerprint,
    envir = benchmark_semantic_fingerprint_cache
  )
  fingerprint
}

benchmark_semantically_equal <- function(
  left,
  right,
  comparison_contract = "strict_graph_isomorphism"
) {
  stopifnot(length(left) == length(right))
  left_key <- benchmark_comparison_key(left, comparison_contract)
  right_key <- benchmark_comparison_key(right, comparison_contract)
  equivalent <- !is.na(left_key) &
    !is.na(right_key) &
    nzchar(left_key) &
    nzchar(right_key) &
    left_key == right_key
  candidates <- which(
    !equivalent &
      !is.na(left) &
      !is.na(right) &
      nzchar(left) &
      nzchar(right)
  )
  if (length(candidates) == 0L) {
    return(equivalent)
  }
  inventory_matches <- vapply(
    candidates,
    function(index) {
      identical(
        benchmark_semantic_inventory_key(
          left[[index]],
          comparison_contract
        ),
        benchmark_semantic_inventory_key(
          right[[index]],
          comparison_contract
        )
      )
    },
    logical(1)
  )
  candidates <- candidates[inventory_matches]
  if (length(candidates) == 0L) {
    return(equivalent)
  }
  equivalent[candidates] <- vapply(
    candidates,
    function(index) {
      left_fingerprint <- benchmark_semantic_fingerprint(
        left[[index]],
        comparison_contract
      )
      right_fingerprint <- benchmark_semantic_fingerprint(
        right[[index]],
        comparison_contract
      )
      !is.na(left_fingerprint) &&
        !is.na(right_fingerprint) &&
        identical(left_fingerprint, right_fingerprint)
    },
    logical(1)
  )
  equivalent
}

benchmark_tool_version <- function(command, arguments = character()) {
  output <- tryCatch(
    system2(command, arguments, stdout = TRUE, stderr = TRUE),
    error = function(error) conditionMessage(error)
  )
  paste(output, collapse = "\n")
}
