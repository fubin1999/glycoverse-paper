library(glycoverse)
library(tidyverse)

# Load and preprocess the data-----
## Expression matrix-----
abundance <- read_csv("data/benchmark/bao2021_hmo_abundance.csv") |>
  rename(sample = `...1`) |>
  mutate(sample = as.character(sample))
expr_mat <- abundance |>
  column_to_rownames("sample") |>
  as.matrix() |>
  t()

## Sample information-----
sample_info <- read_csv("data/benchmark/bao2021_hmo_annotation.csv") |>
  rename(sample = Dataset, patient = `Pateint ID`, dpp = DPP) |>
  mutate(sample = as.character(sample)) |>
  mutate(group = if_else(patient %in% c("L3", "L2"), "non-secretor", "secretor"))

## Variable information-----
glycoct_files <- fs::dir_ls("data/benchmark/bao2021_hmo_glycoct", glob = "*.glycoct_condensed")
glycoct_strings <- map_chr(glycoct_files, read_file)
glycan_structures <- parse_glycoct(unname(glycoct_strings))
glycan_names <- names(glycoct_strings) |>
  str_split_i(fixed("/"), -1) |>
  str_split_i(fixed("."), 1)
var_info <- tibble(
  variable = glycan_names,
  glycan_composition = as_glycan_composition(glycan_structures),
  glycan_structure = glycan_structures
)

## Pack into an experiment object-----
exp <- experiment(
  expr_mat = expr_mat,
  sample_info = sample_info,
  var_info = var_info,
  exp_type = "glycomics",
  glycan_type = "O-Glc"
)
clean_exp <- exp |> normalize_total_area()

# Motif quantification-----
motifs <- c(
  "X40" = "Gal(b1-3)GlcNAc(b1-3)Gal(b1-4)Glc(?1-",
  "X62" = "Gal(b1-3)[Neu5Ac(a2-6)]GlcNAc(b1-3)Gal(b1-4)Glc(?1-",
  "X106" = "Neu5Ac(a2-3)Gal(b1-3)[Neu5Ac(a2-6)]GlcNAc(b1-3)Gal(b1-4)Glc(?1-",
  "X141" = "Neu5Ac(a2-3)Gal(b1-3)[Neu5Ac(a2-6)]GlcNAc(b1-3)[Gal(b1-4)GlcNAc(b1-6)]Gal(b1-4)Glc(?1-"
)
motif_exp <- quantify_motifs(clean_exp, motifs, alignments = "core")
motif_data <- as_tibble(motif_exp) |> select(-sample, -variable)

ratio_data <- motif_data |>
  select(-motif_structure) |>
  pivot_wider(names_from = motif, values_from = value) |>
  mutate(
    ratio_x62_x40 = X62 / X40,
    ratio_x106_x62 = X106 / X62,
  ) |>
  select(patient, dpp, group, starts_with("ratio"))

motif_plots <- motif_data |>
  mutate(motif = factor(motif, levels = c("X40", "X62", "X106", "X141"))) |>
  ggplot(aes(x = log(dpp), y = value)) +
  geom_point(aes(color = group, shape = patient)) +
  geom_smooth(aes(color = group, fill = group), method = "lm", alpha = 0.1) +
  facet_wrap(~motif, scales = "free_y", nrow = 1) +
  guides(shape = guide_legend(nrow = 1)) +
  labs(x = "log(DPP)", y = "Relative Abundance") +
  scale_fill_manual(values = c("secretor" = "#E26A4D", "non-secretor" = "#7B8AA6")) +
  scale_color_manual(values = c("secretor" = "#E26A4D", "non-secretor" = "#7B8AA6")) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    legend.position = "bottom"
  )
ggimage::ggpreview(plot = motif_plots, width = 8, height = 3)
ggsave("results/figures/bao2021_motif_quant.pdf", motif_plots, width = 8, height = 3)

export_cartoons(motifs, "results/figures/bao2021_motifs")

glycan_plots <- as_tibble(clean_exp) |>
  filter(variable %in% c("LNT", "LSTb", "DSLNT", "DSLNH")) |>
  mutate(variable = factor(variable, levels = c("LNT", "LSTb", "DSLNT", "DSLNH"))) |>
  ggplot(aes(x = log(dpp), y = value)) +
  geom_point(aes(color = group, shape = patient)) +
  geom_smooth(aes(color = group, fill = group), method = "lm", alpha = 0.1) +
  facet_wrap(~variable, scales = "free_y", nrow = 1) +
  guides(shape = guide_legend(nrow = 1)) +
  labs(x = "log(DPP)", y = "Relative Abundance") +
  scale_fill_manual(values = c("secretor" = "#E26A4D", "non-secretor" = "#7B8AA6")) +
  scale_color_manual(values = c("secretor" = "#E26A4D", "non-secretor" = "#7B8AA6")) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    legend.position = "bottom"
  )
ggsave("results/figures/bao2021_glycan_quant.pdf", glycan_plots, width = 8, height = 3)
