library(tidyverse)
library(glycoverse)

set.seed(1)

# Load and prepare the data-----
data <- read_csv("data/benchmark/bcc2018_data.csv") |>
  mutate(
    iupac = str_replace_all(glycan, fixed("OS"), "?S"),
    glycan = paste0("G", row_number()),
    .after = 1
  )

expr_mat <- data |>
  select(-iupac) |>
  column_to_rownames("glycan")

var_info <- data |>
  select(glycan, iupac) |>
  mutate(
    iupac = paste0(iupac, "(a1-"),
    glycan_structure = parse_iupac_condensed(iupac),
    glycan_composition = as_glycan_composition(glycan_structure)
  ) |>
  select(variable = glycan, glycan_composition, glycan_structure)

sample_info <- tibble(sample = colnames(expr_mat)) |>
  mutate(group = factor(str_split_i(sample, "_", 1)))

exp <- experiment(
  expr_mat = expr_mat,
  var_info = var_info,
  sample_info = sample_info,
  exp_type = "glycomics",
  glycan_type = "O-GalNAc"
)
clean_exp <- auto_clean(exp)

# Motif quantification and DEA-----
motif_exp <- quantify_motifs(clean_exp, dynamic_motifs(max_size = 5))

dea_res <- motif_exp |>
  gly_ttest() |>
  get_tidy_result()

dea_res_coda <- motif_exp |>
  auto_coda() |>
  gly_ttest() |>
  get_tidy_result()

set_test_res <- gly_set_test(
  motif_exp
)

set_dea_res <- set_test_res |>
  get_tidy_result("sets") |>
  mutate(
    variables = map_chr(variables, \(x) str_c(x, collapse = ";")),
    estimate = map_chr(estimate, \(x) {
      str_c(names(x), signif(x, 6), sep = "=", collapse = ";")
    })
  )

set_dea_res |>
  filter(p_adj < 0.05)

# Save the results-----
motif_data <- as_tibble(motif_exp)
write_csv(motif_data, "results/data/bcc2018_motif_quant.csv")
write_csv(dea_res, "results/data/bcc2018_motif_dea_results.csv")
write_csv(dea_res_coda, "results/data/bcc2018_motif_coda_dea_results.csv")
write_csv(set_dea_res, "results/data/bcc2018_motif_set_dea_results.csv")

# Plots-----
dea_p <- ggplot(dea_res, aes(effect_size, -log10(p_adj))) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "Motif DEA Results",
    y = "-log10 Adjusted P-value",
    x = "Effect Size"
  )

coda_dea_p <- ggplot(dea_res_coda, aes(effect_size, -log10(p_adj))) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "Motif DEA Results",
    y = "-log10 Adjusted P-value",
    x = "Effect Size"
  )

set_dea_p <- ggplot(set_dea_res, aes(effect_size, -log10(p_adj))) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "Motif DEA Results",
    y = "-log10 Adjusted P-value",
    x = "Effect Size"
  )

ggsave(dea_p, filename = "results/figures/glydet_bcc_dea.pdf", width = 4, height = 4)
ggsave(coda_dea_p, filename = "results/figures/glydet_bcc_coda_dea.pdf", width = 4, height = 4)
ggsave(set_dea_p, filename = "results/figures/glydet_bcc_set_dea.pdf", width = 4, height = 4)

export_cartoons(motif_exp$var_info$motif_structure, "results/figures/glydet_bcc_cartoons")
