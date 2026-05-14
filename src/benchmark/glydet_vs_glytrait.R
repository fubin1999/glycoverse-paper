library(tidyverse)
library(glycoverse)
library(glyanno)
library(glydb)

abundance <- read_csv("data/benchmark/glytrait_abundance.csv")
groups <- read_csv("data/benchmark/glytrait_groups.csv")

expr_mat <- abundance |>
  column_to_rownames("Sample") |>
  as.matrix() |>
  t()

sample_info <- groups |>
  rename(sample = Sample, group = Group) |>
  mutate(group = factor(group))

var_info <- tibble(variable = rownames(expr_mat)) |>
  mutate(
    glycan_composition = as_glycan_composition(variable),
    glycan_structure = comp_to_struc(
      glycan_composition,
      db = glydb_structures(
        structure_level = "topological",
        species = "Homo sapiens",
        glycan_type = "N"
      ),
      return_best = TRUE
    )
  )

exp <- experiment(
  expr_mat,
  sample_info = sample_info,
  var_info = var_info,
  exp_type = "glycomics",
  glycan_type = "N"
)

clean_exp <- exp |>
  filter_var(!is.na(glycan_structure)) |>
  auto_clean()

trait_exp <- derive_traits(clean_exp)
trait_dea_res <- gly_ttest(trait_exp)

trait_dea_res |>
  get_tidy_result() |>
  write_csv("results/data/glydet_trait_dea_results.csv")

boxplots <- trait_exp |>
  filter_sig_vars(trait_dea_res) |>
  plot_boxplot()
ggsave("results/figures/glydet_trait_boxplots.pdf", plot = boxplots, width = 6, height = 6)
