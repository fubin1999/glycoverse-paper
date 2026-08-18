library(glycoverse)
library(tidyverse)
library(rstatix)

raw_data <- read_csv("data/benchmark/CLRALR_raw.csv") |>
  mutate(`glycan (m/z)` = as.character(`glycan (m/z)`))
expr_mat <- raw_data |> column_to_rownames("glycan (m/z)")
sample_info <- tibble(sample = colnames(expr_mat))
var_info <- tibble(
  variable = rownames(expr_mat),
  glycan_composition = as_glycan_composition(c(
    "Man(9)GlcNAc(2)",
    "Gal(2)Man(3)GlcNAc(4)Fuc(1)",
    "Gal(2)Man(3)GlcNAc(3)Neu5Ac(2)",
    "Gal(2)Man(3)GlcNAc(4)Neu5Ac(2)",
    "Gal(2)Man(3)GlcNAc(3)Neu5Ac(2)Fuc(1)"
  ))
)
exp <- experiment(expr_mat, sample_info, var_info, exp_type = "glycomics", glycan_type = "N")

rela_exp <- exp |>
  impute_miss_forest() |>
  normalize_total_area()

clean_exp <- auto_clean(exp)

glycowork_data <- read_csv("data/benchmark/CLRALR_glycowork.csv") |>
  mutate(glycan_composition = as_glycan_composition(c(
    "Man(9)GlcNAc(2)",
    "Gal(2)Man(3)GlcNAc(4)Fuc(1)",
    "Gal(2)Man(3)GlcNAc(3)Neu5Ac(2)",
    "Gal(2)Man(3)GlcNAc(4)Neu5Ac(2)",
    "Gal(2)Man(3)GlcNAc(3)Neu5Ac(2)Fuc(1)"
  ))) |>
  mutate(variable = as.character(`glycan (m/z)`), .keep = "unused") |>
  pivot_longer(-c(glycan_composition, variable), names_to = "sample", values_to = "value")

plot_data <- bind_rows(
  as_tibble(exp) |> mutate(type = "raw"),
  as_tibble(rela_exp) |> mutate(type = "relative"),
  as_tibble(clean_exp) |> mutate(type = "glyclean"),
  glycowork_data |> mutate(type = "glycowork"),
) |>
  mutate(group = str_split_i(sample, "_", 1L)) |>
  mutate(
    group = factor(group, levels = c("2pmol", "6pmol", "12pmol", "18pmol")),
    type = factor(type, levels = c("raw", "relative", "glycowork", "glyclean"))
  )

bar_plot <- ggplot(plot_data, aes(as.character(glycan_composition), value)) +
  geom_bar(aes(fill = group), stat = "summary", fun = "mean", position = position_dodge(width = 0.8), width = 0.8) +
  geom_point(aes(group = group), position = position_dodge(width = 0.8), alpha = 0.3) +
  scale_fill_brewer(palette = "Oranges") +
  facet_wrap(~ type, ncol = 1, scale = "free_y") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank(),
  )
ggimage::ggpreview(plot = bar_plot, width = 6, height = 6)
ggsave("results/figures/glyclean_barplot.pdf", bar_plot, width = 6, height = 6)
