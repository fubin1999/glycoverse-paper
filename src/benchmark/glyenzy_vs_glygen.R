library(tidyverse)
library(glycoverse)
library(glydb)
library(glyenzy)

glygen_data <- read_csv("data/benchmark/glygen_glycan_enzyme.csv") |>
  filter(species == "Homo sapiens") |>
  select(glytoucan_ac, gene_name)

glygen_enzyme_genes <- unique(glygen_data$gene_name)

glygen_glycans <- unique(glygen_data$glytoucan_ac)
glydb_glycans <- glydb::glydb_data |>
  filter(
    glycan_type %in% c("N", "O-GalNAc", "O-Man", "O-Fuc", "O-Glc"),
    str_detect(species, fixed("Homo sapiens")),
    remove_substituents(glycan_structure) == glycan_structure,
    count_mono(glycan_composition, "Neu5Gc") == 0
  ) |>
  pull(glytoucan_ac)
common_glycans <- intersect(glygen_glycans, glydb_glycans)

structures <- glydb_data |>
  filter(glytoucan_ac %in% common_glycans) |>
  pull(glycan_structure)

glyenzy_anno <- find_enzyme(structures)
names(glyenzy_anno) <- glydb_data |>
  filter(glytoucan_ac %in% common_glycans) |>
  pull(glytoucan_ac)

glyenzy_anno_df <- glyenzy_anno |>
  enframe(name = "glytoucan_ac", value = "enzyme_gene") |>
  unnest(enzyme_gene) |>
  mutate(glyenzy = TRUE) |>
  pivot_wider(names_from = glytoucan_ac, values_from = glyenzy) |>
  pivot_longer(-enzyme_gene, names_to = "glytoucan_ac", values_to = "glyenzy") |>
  mutate(glyenzy = if_else(is.na(glyenzy), FALSE, TRUE))

glygen_anno_df <- glygen_data |>
  filter(glytoucan_ac %in% common_glycans) |>
  distinct() |>
  rename(enzyme_gene = gene_name) |>
  mutate(glygen = TRUE) |>
  pivot_wider(names_from = glytoucan_ac, values_from = glygen) |>
  pivot_longer(-enzyme_gene, names_to = "glytoucan_ac", values_to = "glygen") |>
  mutate(glygen = if_else(is.na(glygen), FALSE, TRUE))

anno_df <- inner_join(glyenzy_anno_df, glygen_anno_df, by = c("enzyme_gene", "glytoucan_ac"))

write_csv(anno_df, "results/data/glyenzy_vs_glygen.csv")

anno_df |>
  summarise(accuracy = mean(glyenzy == glygen))

conf_mat <- anno_df |>
  group_by(glyenzy, glygen) |>
  summarise(count = n(), .groups = "drop") |>
  ggplot(aes(x = glyenzy, y = glygen, fill = count)) +
  geom_tile() +
  geom_text(aes(label = count), color = "black") +
  scale_fill_distiller("Blues", direction = 1) +
  coord_equal() +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid = element_blank()
  )
ggsave("results/figures/glyenzy_vs_glygen_conf_mat.pdf", conf_mat, width = 3, height = 3)
