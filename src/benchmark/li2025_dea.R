library(tidyverse)
library(glycoverse)

# Load and prepare data-----
data <- read_csv("data/benchmark/li2025_data.csv")

expr_mat <- data |>
  pivot_wider(names_from = "sample", values_from = "value") |>
  column_to_rownames("glycan") |>
  as.matrix()

sample_info <- tibble(sample = colnames(expr_mat)) |>
  mutate(group = sample |> str_remove_all("\\d+") |> as.factor())

var_info <- tibble(variable = rownames(expr_mat)) |>
  mutate(glycan_composition = as_glycan_composition(variable))

exp <- experiment(
  expr_mat = expr_mat,
  sample_info = sample_info,
  var_info = var_info,
  exp_type = "glycomics",
  glycan_type = "N"
) |>
  filter_obs(group != "QC")

clean_exp <- auto_clean(exp)
clean_exp |>
  get_expr_mat() |>
  as.data.frame() |>
  rownames_to_column("glycan") |>
  write_csv("results/data/li2025_clean_data.csv")

# DEA-----
dea_res <- gly_anova(clean_exp)
coda_dea_res <- clean_exp |>
  auto_coda() |>
  gly_anova()

dea_res |>
  get_tidy_result("main_test") |>
  write_csv("results/data/li2025_main_test.csv")

coda_dea_res |>
  get_tidy_result("main_test") |>
  write_csv("results/data/li2025_coda_main_test.csv")

dea_res |>
  get_tidy_result("post_hoc_test") |>
  write_csv("results/data/li2025_post_hoc.csv")

coda_dea_res |>
  get_tidy_result("post_hoc_test") |>
  write_csv("results/data/li2025_coda_post_hoc.csv")

# Compare to orginal results-----
original_res <- read_csv("data/benchmark/li2025_dea_res.csv") |>
  filter(feature_type == "glycan")

original_diff_glycans <- original_res |>
  filter(p.adj < 0.05) |>
  pull(feature)

dea_diff_glycans <- dea_res |>
  get_tidy_result("main_test") |>
  filter(p_adj < 0.05) |>
  pull(variable)

length(intersect(original_diff_glycans, dea_diff_glycans))

coda_dea_diff_glycans <- coda_dea_res |>
  get_tidy_result("main_test") |>
  filter(p_adj < 0.05) |>
  pull(variable)

length(intersect(dea_diff_glycans, coda_dea_diff_glycans))

# Compare to glycowork results-----
glycowork_res <- read_csv("results/data/li2025_glycowork_main_test.csv")

glycowork_diff_glycans <- glycowork_res |>
  filter(significant) |>
  pull(Glycan)

length(intersect(dea_diff_glycans, glycowork_diff_glycans))
