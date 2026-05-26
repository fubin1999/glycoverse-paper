library(tidyverse)
library(glycoverse)
library(glyanno)
library(glydb)
library(ggsignif)

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

plot_data <- trait_exp |>
  filter_var(trait %in% p_df$trait) |>
  as_tibble()

p_df <- trait_dea_res |>
  get_tidy_result() |>
  filter(p_adj < 0.05) |>
  select(trait, p_adj) |>
  mutate(start = "Control", end = "HCC") |>
  mutate(p_adj = paste0("p = ", scales::label_scientific()(p_adj)))

y_pos_df <- plot_data |>
  summarise(y_pos = max(value) + 0.1 * (max(value) - min(value)), .by = trait)

anno_df <- p_df |>
  left_join(y_pos_df, by = "trait")

boxplots <- ggplot(plot_data, aes(group, value)) +
  geom_boxplot(aes(color = group)) +
  facet_wrap(~trait, scales = "free_y") +
  geom_signif(
    data = anno_df,
    mapping = aes(xmin = start, xmax = end, annotations = p_adj, y_position = y_pos),
    manual = TRUE,
    vjust = -0.2,
    size = 0.3,
    textsize = 3.5
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) +
  labs(x = "Group", y = "Trait Value") +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    panel.grid = element_blank()
  )

ggimage::ggpreview(plot = boxplots, width = 6, height = 5)
ggsave("results/figures/glydet_trait_boxplots.pdf", plot = boxplots, width = 6, height = 5)
