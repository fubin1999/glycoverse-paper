library(tidyverse)
library(glycoverse)
library(ggimage)
library(patchwork)

# Read data and create experiment----
expr_mat <- read.csv("data/glycomics/abundance.csv", row.names = 1) |>
  as.matrix() |> t()
var_info <- read_csv("data/glycomics/glycan_structures.csv") |>
  mutate(
    glycan_structure = parse_glycoct(structure),
    glycan_composition = as_glycan_composition(glycan_structure)
  ) |>
  select(variable = composition, glycan_composition, glycan_structure)
sample_info <- read_csv("data/glycomics/groups.csv") |>
  filter(group != "QC") |>
  mutate(group = factor(group, levels = c("HC", "CHB", "LC", "HCC")))

common_samples <- intersect(colnames(expr_mat), sample_info$sample)
common_variables <- intersect(rownames(expr_mat), var_info$variable)

expr_mat <- expr_mat[common_variables, common_samples]
var_info <- var_info |>
  filter(variable %in% common_variables)
sample_info <- sample_info |>
  filter(sample %in% common_samples)

exp <- experiment(
  expr_mat = expr_mat,
  var_info = var_info,
  sample_info = sample_info,
  exp_type = "glycomics",
  glycan_type = "N"
)
set.seed(123)
clean_exp <- auto_clean(exp)

write_rds(clean_exp, "results/data/glycomics_clean_exp.rds")
clean_exp <- read_rds("results/data/glycomics_clean_exp.rds")

# Derived trait analysis----
trait_exp <- derive_traits(clean_exp)
trait_kruskal_res <- gly_kruskal(trait_exp)

trait_exp |>
  filter_sig_vars(trait_kruskal_res) |>
  as_tibble() |>
  mutate(trait = case_match(trait,
    "CA2" ~ "Bi-Antennary",
    "CA3" ~ "Tri-Antennary",
    "CA4" ~ "Tetra-Antennary",
    "GA" ~ "Gal. per Ant.",
    "SG" ~ "Sia. per Ant.",
    "TB" ~ "Bisecting",
    "TC" ~ "Complex",
    "TF" ~ "Fucosylation",
    "TFa" ~ "Arm-Fucosylation",
    "TFc" ~ "Core-Fucosylation",
    "TH" ~ "Hybrid",
    "TS" ~ "Sialylation",
  )) |>
  ggplot(aes(group, value)) +
  geom_boxplot(aes(fill = group, color = group), alpha = 0.5) +
  scale_fill_manual(values = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07")) +
  scale_color_manual(values = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07")) +
  facet_wrap(~trait, nrow = 3, scales = "free_y") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    strip.background = element_blank(),
    axis.ticks.x = element_blank(),
  )
ggsave("results/figures/glycomics_trait_boxplots.pdf", width = 6, height = 5)

# Motif analysis----
motifs <- c(
  GlcNAc = "GlcNAc(??-?)Man(??-",
  LacNAc = "Gal(??-?)GlcNAc(??-?)Man(??-",
  `Sialyl-LacNAc` = "Neu5Ac(??-?)Gal(??-?)GlcNAc(??-?)Man(??-"
)

motif_exp <- quantify_motifs(clean_exp, motifs)

motif_data <- motif_exp |>
  as_tibble() |>
  select(-variable)

# Plot Motifs with individual boxplots stitched together----
# Helper to draw a single motif boxplot
plot_motif_box <- function(motif_name) {
  motif_data |>
    filter(motif == motif_name) |>
    ggplot(aes(group, value, color = group, fill = group)) +
    geom_boxplot(alpha = 0.5) +
    scale_color_manual(values = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07")) +
    scale_fill_manual(values = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07")) +
    labs(
      x = NULL,
      y = "Motif Quantification",
      title = motif_name,
      fill = NULL,
      color = NULL
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.position = "none",
      strip.background = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
}

motif_plots <- purrr::map(unique(motif_data$motif), plot_motif_box)

combined_motif_plot <- wrap_plots(motif_plots, nrow = 1) +
  plot_annotation(tag_levels = "a")

ggpreview(plot = combined_motif_plot, width = 7.5, height = 2.8)
ggsave("results/figures/glycomics_motifs.pdf", combined_motif_plot, width = 7.5, height = 2.8)

# Calculate Flux----
flux_data <- motif_exp |>
  as_tibble() |>
  select(-variable) |>
  pivot_wider(names_from = motif, values_from = value) |>
  mutate(
    flux_lacnac = LacNAc / GlcNAc,
    flux_sia_lacnac = `Sialyl-LacNAc` / LacNAc,
    .keep = "unused"
  ) |>
  pivot_longer(cols = -c(sample, group), names_to = "flux", values_to = "value", names_prefix = "flux_")

kruskal_res <- flux_data |>
  group_by(flux) |>
  kruskal_test(value ~ group)

dunn_res <- flux_data |>
  group_by(flux) |>
  dunn_test(value ~ group)

# Plot Flux----
flux_data |>
  filter(flux %in% c("lacnac", "sia_lacnac")) |>
  mutate(group_int = as.integer(group)) |>
  ggplot(aes(group_int, value, color = flux, fill = flux)) +
  geom_smooth(alpha = 0.1) +
  scale_x_continuous(labels = c("HC", "CHB", "LC", "HCC"), name = NULL) +
  scale_color_manual(values = c("#F2BE30", "#9B5683")) +
  scale_fill_manual(values = c("#F2BE30", "#9B5683")) +
  labs(y = "Flux", fill = NULL, color = NULL) +
  theme_classic()
ggpreview(width = 3.5, height = 2)

ggsave("results/figures/glycomics_flux.pdf", width = 3.5, height = 2)
