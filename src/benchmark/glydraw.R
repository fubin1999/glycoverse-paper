library(glydraw)
library(glyparse)
library(glyrepr)
library(tidyverse)

seqs <- read_csv("data/benchmark/text_nomenclatures/glycan_sequences_wurcs.csv")

set.seed(123)
plot_data <- seqs |>
  slice_sample(n = 2000) |>
  mutate(glycan = parse_wurcs(sequence_wurcs, on_failure = "na", progress = TRUE, drop_generic = TRUE)) |>
  filter(!is.na(glycan)) |>
  slice_head(n = 1000)

plots <- plot_data |>
  mutate(
    fp = fs::path("results/figures/glydraw_cartoons/", glytoucan_ac, ext = "png"),
    cartoon = map(glycan, draw_cartoon, .progress = TRUE)
  )
walk2(plots$cartoon, plots$fp, save_cartoon, .progress = TRUE)
