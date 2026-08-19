#!/usr/bin/env Rscript

# Reproducible entry point for the glyparse corpus benchmark.
#
# Usage:
#   Rscript src/benchmark/glyparse/glyparse.R
#   Rscript src/benchmark/glyparse/glyparse.R --stage glyparse
#   Rscript src/benchmark/glyparse/glyparse.R --stage external
#   Rscript src/benchmark/glyparse/glyparse.R --stage adjudicate
#   Rscript src/benchmark/glyparse/glyparse.R --stage report

script_path <- sub(
  "^--file=",
  "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1L]]
)
script_dir <- dirname(normalizePath(script_path, mustWork = TRUE))
module_dir <- script_dir

arguments <- commandArgs(trailingOnly = TRUE)
stage <- "all"
if (length(arguments) > 0L) {
  if (!identical(arguments[[1L]], "--stage") || length(arguments) != 2L) {
    stop(
      paste0(
        "Usage: Rscript src/benchmark/glyparse/glyparse.R ",
        "[--stage glyparse|external|adjudicate|report]"
      )
    )
  }
  stage <- arguments[[2L]]
}
allowed <- c("all", "glyparse", "external", "adjudicate", "report")
if (!stage %in% allowed) {
  stop("Unknown stage: ", stage)
}

run_r <- function(name) {
  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    file.path(module_dir, name)
  )
  if (!identical(status, 0L)) {
    stop(name, " failed with exit status ", status)
  }
}

if (stage %in% c("all", "glyparse")) {
  run_r("run-glyparse.R")
}
if (stage %in% c("all", "external")) {
  run_r("run-external.R")
}
if (stage %in% c("all", "adjudicate")) {
  run_r("adjudicate.R")
}
if (stage %in% c("all", "report")) {
  run_r("classify-disagreements.R")
  run_r("build-report.R")
  run_r("build-disagreement-report.R")
}
