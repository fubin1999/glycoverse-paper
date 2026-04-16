library(glycoverse)
library(tidyverse)
library(glydb)

motif_data <- read_csv("data/benchmark/glycan_motif.csv")
glycans <- glydb_data$glycan_structure
glymotif_result <- have_motifs(glycans, db_motifs(), strict_sub = FALSE)

# Export glydb data for glycowork
write_csv(glydb_data, "data/benchmark/glydb_data.csv")
write_csv(glymotif:::glygen_motifs, "data/benchmark/glygen_motifs.csv")

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

write_csv(res_df, "results/data/glymotif_vs_glygen.csv")

sens_data <- res_df |>
  filter(glygen) |>
  summarise(sensitivity = mean(glymotif), .by = motif_name)

sens_barplot <- ggplot(sens_data, aes(motif_name, sensitivity)) +
  geom_col(fill = "orange") +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = "Motif Name",
    y = "Sensitivity",
    title = "Glymotif Sensitivity for GlyGen Motifs"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.ticks.x = element_blank(),
  )
ggimage::ggpreview(plot = sens_barplot, width = 12, height = 4)
ggsave("results/figures/glymotif_sensitivity_barplot.pdf", sens_barplot, width = 12, height = 4)

res_df |>
  filter(glygen, !glymotif)
# One false positive: G45714BQ doesn't have motif GD1a
