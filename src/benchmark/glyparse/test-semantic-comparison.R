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

message("Semantic comparison regression checks passed.")
