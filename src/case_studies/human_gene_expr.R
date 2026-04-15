library(tidyverse)
library(edgeR)
library(limma)
library(ComplexHeatmap)
library(ggimage)
library(cowplot)

data <- read_tsv("data/glycomics/GSE210545_gene_expression.txt")

bubble_data <- data |>
  filter(str_starts(gene_symbol, "B4GALT") | str_detect(gene_symbol, "ST[36]GAL\\d+")) |>
  select(gene_symbol, starts_with("fpkm")) |>
  pivot_longer(-gene_symbol, names_to = "sample", values_to = "fpkm") |>
  mutate(
    group = str_split_i(sample, "_", 2),
    group = factor(group, levels = c("Normal", "Disease")),
    sample_id = str_split_i(sample, "group", 2),
    sample = paste0(group, sample_id),
  ) |>
  mutate(
    log_fpkm = log2(fpkm + 1),
    z_score = as.numeric(scale(log_fpkm)),
    .by = gene_symbol
  ) |>
  summarise(mean_z_score = mean(z_score), .by = c(gene_symbol, group)) |>
  mutate(
    enzyme_type = if_else(str_starts(gene_symbol, "B4GALT"), "B4GALT", "ST3/6GAL"),
    enzyme_type = factor(enzyme_type),
  )

# differential analysis for B4GALT and ST3/6GAL genes using limma
glyco_sample_cols <- names(data) |> keep(~str_starts(.x, "read_count_"))
glyco_samples <- str_remove(glyco_sample_cols, "^read_count_")
glyco_groups <- str_split_i(glyco_samples, "_", 1) |>
  factor(levels = c("Normal", "Disease"))

limma_matrix <- data |>
  filter(str_starts(gene_symbol, "B4GALT") | str_detect(gene_symbol, "ST[36]GAL\\d+")) |>
  select(gene_symbol, all_of(glyco_sample_cols)) |>
  mutate(across(all_of(glyco_sample_cols), ~log2(.x + 1))) |>
  column_to_rownames("gene_symbol") |>
  as.matrix()

design <- model.matrix(~0 + glyco_groups)
colnames(design) <- c("Normal", "Disease")

limma_fit <- lmFit(limma_matrix, design)
limma_contrast <- makeContrasts(DiseaseVsNormal = Disease - Normal, levels = design)
limma_fit2 <- limma_fit |>
  contrasts.fit(limma_contrast) |>
  eBayes()

limma_results <- topTable(limma_fit2, number = Inf, sort.by = "none") |>
  rownames_to_column("gene_symbol") |>
  mutate(
    enzyme_type = if_else(str_starts(gene_symbol, "B4GALT"), "B4GALT", "ST3/6GAL")
  )

write_tsv(limma_results, "results/data/glycomics_gene_limma.tsv")

plot_bubble <- function(data) {
  ggplot(data, aes(group, fct_rev(gene_symbol), size = mean_z_score)) +
    geom_point(aes(color = group)) +
    scale_color_manual(values = c(Normal = "#2E9FDF", Disease = "#E7B800")) +
    scale_x_discrete(labels = c("HC", "LC")) +
    facet_wrap(~enzyme_type, nrow = 1) +
    labs(size = "Z-Score") +
    guides(color = "none") +
    theme_minimal() +
    theme(
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
    )
}

stgal_bubble <- bubble_data |>
  filter(enzyme_type == "ST3/6GAL") |>
  plot_bubble()
b4galt_bubble <- bubble_data |>
  filter(enzyme_type == "B4GALT") |>
  plot_bubble()

plot_grid(stgal_bubble, b4galt_bubble, nrow = 1)

ggpreview(width = 4.5, height = 2)
ggsave("results/figures/glycomics_gene_bubble.pdf", width = 4.5, height = 2)

