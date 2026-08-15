#!/usr/bin/env Rscript

script_path <- sub(
  "^--file=",
  "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1L]]
)
source(file.path(dirname(normalizePath(script_path)), "common.R"))

benchmark_load_packages()

glyparse_g12345 <- paste0(
  "Man(?1-?)[Man(??-?)]Man(?1-?)GlcNAc(?1-?)",
  "[Fuc(?1-?)]GlcNAc(?1-"
)
reference_g12345 <- paste0(
  "Man(??-?)[Man(?1-?)]Man(?1-?)GlcNAc(?1-?)",
  "[Fuc(?1-?)]GlcNAc(?1-"
)

stopifnot(
  glyparse_g12345 != reference_g12345,
  benchmark_semantically_equal(glyparse_g12345, reference_g12345),
  benchmark_semantically_equal(
    "Man(?1-?)[Man(?1-?)]Man(?1-?)GlcNAc(?1-?)[Fuc(?1-?)]GlcNAc(?1-",
    reference_g12345
  ),
  !benchmark_semantically_equal(
    "Gal(b1-3)GalNAc(a1-",
    "Gal(b1-3)GalNAc(?1-"
  ),
  benchmark_semantically_equal(
    "Gal(b1-3)GalNAc(a1-",
    "Gal(b1-3)GalNAc(?1-",
    "graph_isomorphism_ignore_terminal_reducing_end_token"
  ),
  !benchmark_semantically_equal(
    "Gal(b1-3)GalNAc(a1-",
    "Gal(b1-4)GalNAc(a1-"
  )
)

classification_cases <- benchmark_classify_semantic_differences(data.frame(
  difference_reason = c(
    "Residue or residue-substituent composition differs",
    paste(
      "Residue or residue-substituent composition differs",
      "Reducing-end alditol state differs",
      sep = "; "
    ),
    paste(
      "Linkage or anomer multiset differs",
      "Canonical topology or residue-to-residue attachment differs",
      sep = "; "
    )
  ),
  stringsAsFactors = FALSE
))
empty_classification <- benchmark_classify_semantic_differences(data.frame())

stopifnot(
  identical(
    classification_cases$primary_mismatch_class,
    c(
      "residue_composition",
      "reducing_end_alditol",
      "topology_or_attachment"
    )
  ),
  identical(classification_cases$mismatch_class_count, c(1L, 2L, 2L)),
  classification_cases$mismatch_has_residue_composition[[2L]],
  classification_cases$mismatch_has_reducing_end_alditol[[2L]],
  classification_cases$mismatch_has_linkage_or_anomer[[3L]],
  classification_cases$mismatch_has_topology_or_attachment[[3L]],
  nrow(empty_classification) == 0L,
  "primary_mismatch_class" %in% names(empty_classification)
)

message("Semantic comparison regression checks passed.")
