library(glycoverse)
library(tidyverse)

# Prepare data for glymotif-----
# Freeze the cohort to the 6,036 glycans used in the original comparison.
# Do not regenerate this set from the evolving glydb_data object.
glycans <- read_csv(
  "src/benchmark/data/glymotif_glycans.csv",
  col_types = cols(glycan = col_character())
) |>
  pull(glycan)

stopifnot(length(glycans) == 6036L, !anyDuplicated(glycans))

motif_data <- glymotif:::glygen_motifs
motif_structure_column <- intersect(
  c("glycan_structure", "iupac"),
  names(motif_data)
)[[1L]]
motifs <- motif_data[[motif_structure_column]] |>
  as.character() |>
  str_replace("\\(...$", "(??-") |>
  unique()

# Prepare data for glycowork-----
to_glycowork <- function(x) {
  x |>
    str_sub(1L, -5L) |>
    # glycowork uses OS, rather than ?S, for an unknown sulfate position.
    str_replace_all(fixed("?S"), "OS")
}

glycan_map <- tibble(
  glycan = glycans,
  glycowork_glycan = to_glycowork(glycans)
)
motif_map <- tibble(
  motif = motifs,
  glycowork_motif = to_glycowork(motifs)
)

stopifnot(
  identical(
    to_glycowork(c("Gal3S(??-", "Gal3/6S(??-", "Gal?S(??-")),
    c("Gal3S", "Gal3/6S", "GalOS")
  )
)

write_csv(
  transmute(glycan_map, glycan = glycowork_glycan),
  "data/benchmark/glycans_glycowork.csv"
)
write_csv(
  transmute(motif_map, motif = glycowork_motif),
  "data/benchmark/motifs_glycowork.csv"
)

# Run glymotif-----
glymotif_result <- count_motifs(glycans, motifs)
colnames(glymotif_result) <- motifs

# Merge results-----
glycowork_annotations <- read_csv(
  "data/benchmark/glycowork_motif_annotations.csv"
)
names(glycowork_annotations)[[1L]] <- "glycowork_glycan"

glycowork_result <- glycowork_annotations |>
  pivot_longer(
    -glycowork_glycan,
    names_to = "glycowork_motif",
    values_to = "glycowork"
  ) |>
  inner_join(glycan_map, by = join_by(glycowork_glycan)) |>
  inner_join(motif_map, by = join_by(glycowork_motif)) |>
  select(glycan, motif, glycowork)

result_data <- glymotif_result |>
  as.data.frame() |>
  rownames_to_column("glycan") |>
  pivot_longer(-glycan, names_to = "motif", values_to = "glymotif") |>
  inner_join(glycowork_result, by = join_by(glycan, motif))

result_data |>
  filter(glymotif != glycowork)

write_csv(result_data, "results/data/glymotif_vs_glycowork.csv")
