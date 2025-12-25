library(tidyverse)
library(glycoverse)
library(ggradar)
library(ggridges)
library(ggimage)
library(cowplot)

conflicted::conflicts_prefer(stats::aggregate)
conflicted::conflicts_prefer(dplyr::filter)

colors <- c(
  brain = "#2E9FDF",
  heart = "#00BB0C",
  liver = "#E7B800",
  lung = "#FC4E07",
  kidney = "#8B5A2B"
)

# Data import and preprocessing-----
exp <- read_strucgp("data/glycoproteomics/strucgp_result.xlsx")

write_rds(exp, "results/data/strucgp_result.rds")
exp <- read_rds("results/data/strucgp_result.rds")

data <- exp |>
  as_tibble() |>
  filter(value == 1) |>
  select(-value)
sample_info <- read_csv("data/glycoproteomics/sample_info.csv") |>
  mutate(
    sample = str_split_i(Name, fixed("."), 1),
    protease = case_when(
      Trypsin == "+" ~ "Trypsin",
      `Trypsin+Glu-C` == "+" ~ "Trypsin+Glu-C",
      `Trypsin+Lysine-C` == "+" ~ "Trypsin+Lysine-C",
      .default = "others"
    ),
    lc = case_when(
      `Sepharose CL-4B` == "+" ~ "Sepharose CL-4B",
      `ZIC-HILIC` == "+" ~ "ZIC-HILIC",
      .default = "others"
    ),
    organ = case_when(
      Heart == "+" ~ "heart",
      liver == "+" ~ "liver",
      Brain == "+" ~ "brain",
      lung == "+" ~ "lung",
      kidney == "+" ~ "kidney",
      .default = "others"
    ),
    condition = str_c(organ, "_", protease, "_", lc),
    protease = factor(protease),
    lc = factor(lc),
    organ = factor(organ),
    condition = factor(condition),
    .keep = "unused"
  )
data <- data |>
  left_join(sample_info, by = "sample") |>
  distinct(organ, protein, protein_site, gene, glycan_composition, glycan_structure)

# Data import from other sources-----
## pGlyco3-----
pglyco_exp <- read_pglyco3("data/glycoproteomics/pglyco3_result.txt")
pglyco_ident <- pglyco_exp |>
  get_var_info() |>
  distinct(protein, protein_site, glycan_composition)

## Glyco-Decipher-----
decipher_data <- read_csv("data/glycoproteomics/glyco_decipher_result.csv")
decipher_ident <- decipher_data |>
  distinct(Protein, GlycoSite, GlycanComposition) |>
  rename(protein = Protein, protein_site = GlycoSite, glycan_composition = GlycanComposition) |>
  mutate(
    protein = str_remove(protein, ";$"),
    protein_site = str_remove(protein_site, ";$"),
  ) |>
  separate_longer_delim(c(protein, protein_site), delim = ";") |>
  distinct() |>
  mutate(
    protein = str_split_i(protein, fixed("|"), 2),
    protein_site = as.integer(protein_site),
    glycan_composition = str_replace(glycan_composition, "Fuc", "dHex"),
    glycan_composition = as_glycan_composition(glycan_composition)
  )

## MSFragger-----
msfragger_data1 <- read_tsv("data/glycoproteomics/msfragger_results/glu-c/psm.tsv") |>
  janitor::clean_names()
msfragger_data2 <- read_tsv("data/glycoproteomics/msfragger_results/lys-c/psm.tsv") |>
  janitor::clean_names()
msfragger_data3 <- read_tsv("data/glycoproteomics/msfragger_results/trypsin/psm.tsv") |>
  janitor::clean_names()

msfragger_data <- bind_rows(msfragger_data1, msfragger_data2, msfragger_data3)
rm(msfragger_data1, msfragger_data2, msfragger_data3)

msfragger_ident <- msfragger_data |>
  filter(
    str_count(assigned_modifications, "N") == 1,
    str_detect(observed_modifications, "Hex|HexNAc")
  ) |>
  mutate(protein_site = protein_start + as.integer(str_extract(assigned_modifications, "(\\d+)N", group = 1)) - 1) |>
  distinct(protein_id, protein_site, observed_modifications) |>
  rename(protein = protein_id, glycan_composition = observed_modifications) |>
  mutate(
    glycan_composition = str_split_i(glycan_composition, fixed(" % "), 1),
    glycan_composition = str_replace(glycan_composition, "Fuc", "dHex"),
    glycan_composition = as_glycan_composition(glycan_composition)
  )

# Get common glycoforms-----
data <- data |>
  semi_join(pglyco_ident, by = c("protein", "protein_site", "glycan_composition")) |>
  semi_join(decipher_ident, by = c("protein", "protein_site", "glycan_composition")) |>
  semi_join(msfragger_ident, by = c("protein", "protein_site", "glycan_composition"))

# Branching motif analysis-----
branching_motifs <- c(
  "Hex(??-?)[NeuGc(??-?)]HexNAc(??-?)",
  "NeuGc(??-?)Hex(??-?)[NeuGc(??-?)]HexNAc(??-?)",
  "NeuGc(??-?)Hex(??-?)HexNAc(??-?)",
  "NeuGc(??-?)Hex(??-?)[dHex(??-?)]HexNAc(??-?)",
  "NeuAc(??-?)Hex(??-?)[NeuGc(??-?)]HexNAc(??-?)",
  "NeuAc(??-?)Hex(??-?)HexNAc(??-?)",

  "HexNAc(??-?)[dHex(??-?)]HexNAc(??-?)",
  "dHex(??-?)HexNAc(??-?)[dHex(??-?)]HexNAc(??-?)",
  "Hex(??-?)[dHex(??-?)]HexNAc(??-?)",
  "HexNAc(??-?)",
  "NeuAc(??-?)HexNAc(??-?)HexNAc(??-?)",
  "NeuAc(??-?)HexNAc(??-?)[dHex(??-?)]HexNAc(??-?)",

  "Hex(??-?)HexNAc(??-?)",
  "HexNAc(??-?)HexNAc(??-?)",
  "NeuAc(??-?)NeuAc(??-?)Hex(??-?)[NeuAc(??-?)]HexNAc(??-?)",
  "NeuAc(??-?)Hex(??-?)[NeuAc(??-?)]HexNAc(??-?)",
  "Hex(??-?)[NeuAc(??-?)]HexNAc(??-?)",
  "NeuAc(??-?)NeuAc(??-?)Hex(??-?)HexNAc(??-?)"
)
branching_motifs <- str_c(branching_motifs, "Hex(??-?)Hex(??-?)HexNAc(??-?)HexNAc(??-")
names(branching_motifs) <- str_c("m", seq_along(branching_motifs))

branching_motif_res <- have_motifs(data$glycan_structure, branching_motifs, ignore_linkages = TRUE, alignments = "terminal") |> as_tibble()
data_with_branching_motifs <- bind_cols(data, branching_motif_res)

branching_motif_ggradar_data <- data_with_branching_motifs |>
  pivot_longer(all_of(starts_with("m")), names_to = "motif", values_to = "has_motif") |>
  summarise(prop = mean(has_motif), .by = c(organ, motif)) |>
  mutate(prop = prop * 100) |>
  mutate(prop = scales::rescale(prop), .by = motif) |>
  pivot_wider(names_from = motif, values_from = prop) |>
  rename(group = organ)
ggradar(
  branching_motif_ggradar_data,
  values.radar = character(),
  background.circle.transparency = 0,
  group.point.size = 5,
  group.colours = colors[levels(branching_motif_ggradar_data$group)]
)
ggpreview(width = 10, height = 6)
ggsave("results/figures/strucgp_motifs.pdf", width = 10, height = 6)

# Meta-property analysis-----
my_mps <- all_mp_fns()
my_mps[["nAc"]] <- function(x) count_mono(x, "NeuAc")
my_mps[["nGc"]] <- function(x) count_mono(x, "NeuGc")
my_mps[["nS"]] <- function(x) count_mono(x, "NeuAc") + count_mono(x, "NeuGc")
my_mps[["nGN"]] <- function(x) count_motif(x, "HexNAc(??-?)HexNAc(??-")

mp_res <- get_meta_properties(data$glycan_structure, mp_fns = my_mps)
data_with_mps <- bind_cols(data, mp_res)

mp_summary_data <- data_with_mps |>
  summarise(
    `Core Fucosylation` = mean(nFc),
    `Bisecting` = mean(B),
    `Complex Type` = mean(Tp == "complex"),
    `Hybrid Type` = mean(Tp == "hybrid"),
    `High-Mannose Type` = mean(Tp == "highmannose"),
    `Pauci-Mannose Type` = mean(Tp == "paucimannose"),
    `Bi-Antennary` = mean(nA == 2L),
    `Tri-Antennary` = mean(nA == 3L),
    `Tetra-Antennary` = mean(nA == 4L),
    .by = organ
  ) |>
  pivot_longer(cols = -organ, names_to = "mp", values_to = "value") |>
  mutate(mp = factor(mp, levels = c(
    "Core Fucosylation",
    "Bisecting",
    "Complex Type",
    "Hybrid Type",
    "High-Mannose Type",
    "Pauci-Mannose Type",
    "Bi-Antennary",
    "Tri-Antennary",
    "Tetra-Antennary"
  )))

ggplot(mp_summary_data, aes(x = organ, y = value)) +
  geom_col(fill = "grey90", width = 0.1) +
  geom_point(aes(color = organ, size = value)) +
  facet_wrap(~mp, nrow = 3, scales = "free_y") +
  scale_color_manual(values = colors) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.1))) +
  labs(y = "Proportion") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    strip.background = element_blank(),
  )
ggpreview(width = 6, height = 6)
ggsave("results/figures/strucgp_mps.pdf", width = 6, height = 6)

# Maturation analysis-----
maturation_scores <- data_with_mps |>
  mutate(
    type_score = case_match(
      Tp,
      "complex" ~ 1,
      "hybrid" ~ 0.5,
      "highmannose" ~ 0,
      "paucimannose" ~ 0,
    ),
    antennary_score = case_match(
      nA,
      2L ~ 0,
      3L ~ 0.5,
      4L ~ 1,
      .default = 0
    ),
    core_fuc_score = nFc,
    bisecting_score = if_else(B, 1, 0),
    elongation_score = if_else(nA > 0, (nG + nGN) / nA, 0),
    sia_score = if_else(nA > 0, nS / nA, 0),
    arm_fuc_score = if_else(nA > 0, nFa / nA, 0),
  ) |>
  mutate(score = type_score + antennary_score + core_fuc_score + bisecting_score + elongation_score + sia_score + arm_fuc_score)

maturation_scores |>
  ggplot(aes(score, fct_rev(organ), fill = organ)) +
  geom_density_ridges() +
  scale_fill_manual(values = colors) +
  labs(x = "Maturation Score") +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
  )

maturation_scores_summary <- maturation_scores |>
  summarise(
    prop_low = mean(score == 0),
    prop_medium = mean(score > 0 & score <= 3),
    prop_high = mean(score > 3),
    .by = organ
  ) |>
  pivot_longer(-organ, names_to = "score_range", values_to = "prop", names_prefix = "prop_") |>
  mutate(score_range = factor(score_range, levels = c("low", "medium", "high")))

ggplot(maturation_scores_summary, aes(organ, prop)) +
  geom_col(aes(fill = fct_rev(score_range))) +
  geom_text(
    data = maturation_scores_summary |> filter(score_range == "high"),
    aes(label = scales::percent(prop)),
    y = 1.05, size = 3,
    color = RColorBrewer::brewer.pal(3, "Blues")[3]
  ) +
  geom_text(
    data = maturation_scores_summary |>
      filter(score_range == "medium") |>
      left_join(
        maturation_scores_summary |>
          filter(score_range == "low") |>
          select(organ, y = prop) |>
          mutate(y = y - 0.05),
        by = "organ"
      ),
    aes(label = scales::percent(prop, accuracy = 0.1), y = y),
    size = 3,
    color = RColorBrewer::brewer.pal(3, "Blues")[2]
  ) +
  geom_text(
    data = maturation_scores_summary |> filter(score_range == "low"),
    aes(label = scales::percent(prop, accuracy = 0.1)),
    y = 0.05, size = 3,
    color = "white"
  ) +
  scale_fill_brewer(palette = "Blues", direction = -1) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.05)) +
  labs(y = "Proportion", fill = "Maturation\nDegree") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.line.y = element_line(),
    axis.ticks.y = element_line(color = "black"),
    axis.text.y = element_text(color = "black"),
  )
ggpreview(width = 4, height = 2.5)
ggsave("results/figures/strucgp_maturation.pdf", width = 4, height = 2.5)

kruskal_res <- rstatix::kruskal_test(score ~ organ, data = maturation_scores)
dunn_res <- rstatix::dunn_test(score ~ organ, data = maturation_scores, p.adjust.method = "bonferroni")

high_mature_proteins <- maturation_scores |>
  summarise(score = mean(score), .by = c(organ, protein)) |>
  summarise(n_high = sum(score > 3), .by = protein) |>
  filter(n_high > 2) |>
  pull(protein)

low_mature_proteins <- maturation_scores |>
  summarise(score = mean(score), .by = c(organ, protein)) |>
  summarise(n_low = sum(score == 0), .by = protein) |>
  filter(n_low > 2) |>
  pull(protein)

high_go_res <- clusterProfiler::enrichGO(
  gene = high_mature_genes,
  OrgDb = org.Mm.eg.db::org.Mm.eg.db,
  universe = unique(data$gene),
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH"
)
high_go_res <- clusterProfiler::simplify(high_go_res)
high_go_plot <- enrichplot::dotplot(high_go_res)

low_go_res <- clusterProfiler::enrichGO(
  gene = low_mature_genes,
  OrgDb = org.Mm.eg.db::org.Mm.eg.db,
  universe = unique(data$gene),
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH"
)
low_go_res <- clusterProfiler::simplify(low_go_res)
low_go_plot <- enrichplot::dotplot(low_go_res)

plot_grid(high_go_plot, low_go_plot)
ggpreview(width = 10, height = 5)
ggsave("results/figures/strucgp_go.pdf", width = 10, height = 5)

maturation_scores |>
  summarise(
    mean_score = mean(score, na.rm = TRUE),
    n_glycan = n_distinct(glycan_structure),
    .by = c(organ, protein)
  ) |>
  mutate(
    score_range = max(mean_score) - min(mean_score),
    min_n_glycan = min(n_glycan),
    .by = protein
  ) |>
  filter(min_n_glycan > 10, score_range > 2)

maturation_scores |>
  filter(protein == "P11438") |>
  ggplot(aes(score, fct_rev(organ), fill = organ)) +
  geom_density_ridges(alpha = 0.5) +
  scale_fill_manual(values = colors) +
  scale_y_discrete(expand = expansion(mult = c(0, 0.5))) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.title.y = element_blank(),
  )
ggpreview(width = 3, height = 2.5)
ggsave("results/figures/strucgp_P11438_maturation_density.pdf", width = 3, height = 2.5)

maturation_scores |>
  filter(protein == "P11438") |>
  mutate(protein_site = factor(protein_site)) |>
  summarise(mean_score = mean(score, na.rm = TRUE), .by = c(organ, protein_site)) |>
  mutate(score_color = if_else(mean_score > 3, "white", "black")) |>
  ggplot(aes(protein_site, organ, fill = mean_score)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = scales::number(mean_score, accuracy = 0.1), color = score_color)) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  scale_color_identity() +
  coord_equal() +
  labs(x = "Glycosite", y = "Organ", fill = "Maturation\nScore") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.title.y = element_blank(),
  )
ggpreview(width = 4.5, height = 2.5)
ggsave("results/figures/strucgp_P11438_maturation_heatmap.pdf", width = 4.5, height = 2.5)

data |>
  filter(protein == "P11438", protein_site == 70)

data |>
  filter(protein == "P11438", protein_site == 70) |>
  pull(glycan_structure)
