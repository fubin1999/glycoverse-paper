library(glycoverse)
library(tidyverse)
library(glydb)

# Prepare data for glymotif-----
motif_data <- read_csv("data/benchmark/glycan_motif.csv")
glycans <- glydb_data$glycan_structure |>
  as.character() |>
  str_replace("\\(...$", "(??-") |>
  unique()
motifs <- glymotif:::glygen_motifs |>
  mutate(iupac = str_replace(iupac, "\\(...$", "(??-")) |>
  pull(iupac) |>
  unique()

# Prepare data for glycowork-----
glycans_glycowork <- tibble(glycan = glycans) |>
  mutate(glycan = str_sub(glycan, 1L, -5L))
motifs_glycowork <- tibble(motif = motifs) |>
  mutate(motif = str_sub(motif, 1L, -5L))
write_csv(glycans_glycowork, "data/benchmark/glycans_glycowork.csv")
write_csv(motifs_glycowork, "data/benchmark/motifs_glycowork.csv")

# Run glymotif-----
glymotif_result <- count_motifs(glycans, motifs)
colnames(glymotif_result) <- motifs

# Merge results-----
glycowork_result <- read_csv("data/benchmark/glycowork_motif_annotations.csv") |>
  rename(glycan = ...1) |>
  pivot_longer(-glycan, names_to = "motif", values_to = "glycowork") |>
  mutate(
    glycan = paste0(glycan, "(??-"),
    motif = paste0(motif, "(??-")
  )

result_data <- glymotif_result |>
  as.data.frame() |>
  rownames_to_column("glycan") |>
  pivot_longer(-glycan, names_to = "motif", values_to = "glymotif") |>
  inner_join(glycowork_result, by = join_by(glycan, motif))

result_data |>
  filter(glymotif != glycowork)  # 3,153

write_csv(result_data, "results/data/glymotif_vs_glycowork.csv")
