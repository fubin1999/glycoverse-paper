library(glycoverse)
library(tidyverse)
library(glydb)

motif_data <- read_csv("data/benchmark/glycan_motif.csv")
glycans <- glydb_data$glycan_structure
glymotif_result <- have_motifs(glycans, db_motifs(), strict_sub = FALSE)

res_df <- glymotif_result |>
  as.data.frame() |>
  rownames_to_column("glycan_structure") |>
  as_tibble() |>
  pivot_longer(-glycan_structure, names_to = "motif_name", values_to = "glymotif") |>
  left_join(
    glydb_data |>
      select(glytoucan_ac, glycan_structure) |>
      mutate(glycan_structure = as.character(glycan_structure)),
    join_by(glycan_structure)
  ) |>
  left_join(
    motif_data |>
      select(glytoucan_ac, motif_name) |>
      mutate(glygen = TRUE),
    join_by(motif_name, glytoucan_ac)
  ) |>
  mutate(glygen = if_else(is.na(glygen), FALSE, glygen)) |>
  relocate(glytoucan_ac, glycan_structure, motif_name, glymotif, glygen)

filtered_res_df <- res_df |>
  filter(sum(glygen) > 0, .by = motif_name) |>
  filter(sum(glygen) > 0, .by = glycan_structure)
n_distinct(filtered_res_df$motif_name)  # 107
n_distinct(filtered_res_df$glycan_structure)  #3913

filtered_res_df |>
  filter(glygen, !glymotif)
# One false positive: G45714BQ doesn't have motif GD1a

filtered_res_df |>
  filter(!glygen, glymotif) |>
  select(glytoucan_ac, motif_name)

write_csv(filtered_res_df, "results/data/glymotif_vs_glygen.csv")
