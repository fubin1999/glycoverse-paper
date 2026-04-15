library(tidyverse)
library(glycoverse)
library(glydb)
library(bench)
library(ggimage)

structures <- glydb_structures(glycan_type = "N")

run_derive_traits <- function(n_samples, n_structures) {
  structures <- sample(structures, n_structures, replace = FALSE)
  var_info <- tibble(
    variable = paste0("V", 1:n_structures),
    glycan_composition = as_glycan_composition(structures),
    glycan_structure = structures
  )
  sample_info <- tibble(sample = paste0("S", 1:n_samples))
  expr_mat <- matrix(exp(rnorm(n_samples * n_structures)), nrow = n_structures, ncol = n_samples)
  colnames(expr_mat) <- sample_info$sample
  rownames(expr_mat) <- var_info$variable
  exp <- experiment(
    expr_mat = expr_mat,
    var_info = var_info,
    sample_info = sample_info,
    exp_type = "glycomics",
    glycan_type = "N"
  )
  derive_traits(exp)
}

set.seed(123)
sample_sizes <- seq(100, 500, by = 100)
structure_sizes <- seq(100, 500, by = 100)
results <- bench::press(
  n_samples = sample_sizes,
  n_structures = structure_sizes,
  bench::mark(run_derive_traits(n_samples, n_structures))
)

p <- results |>
  mutate(
    n_structures = as.factor(n_structures),
    n_samples = as.factor(n_samples)
  ) |>
  ggplot(aes(n_structures, median, color = n_samples, group = n_samples)) +
  geom_line() +
  geom_point() +
  labs(x = "Number of Structures", y = "Run Time (s)", color = "Number of Samples") +
  theme_minimal() +
  scale_color_viridis_d()
ggpreview(plot = p, width = 6, height = 4)
ggsave("results/figures/glydet_bench.pdf", width = 6, height = 4)
