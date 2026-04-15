library(tidyverse)
library(ggimage)

data <- read_tsv("data/glycoproteomics/expression_atlas/gene_expr_Cmah.tsv") |>
  rename(dataset = `...1`) |>
  pivot_longer(cols = -dataset, names_to = "organ", values_to = "tpm") |>
  mutate(log2_tpm = log2(tpm + 1)) |>
  mutate(evidence = if_else(
    dataset %in% c("Organism part - Huttlin et al", "Huttlin et al., 2010 - tissues - brain, heart, etc", "Meierhofer et al., 2016 - organism part"),
    "proteomics",
    "transcriptomics",
  )) |>
  mutate(evidence = factor(evidence, levels = c("proteomics", "transcriptomics")))

ggplot(data, aes(x = organ, y = dataset)) +
  geom_point(aes(size = log2_tpm), color = "#2E9FDF") +
  scale_color_distiller(palette = "Spectral") +
  labs(title = "Cmah", size = "log2(TPM + 1)") +
  coord_equal() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    plot.title = element_text(hjust = 0.5),
  )
ggpreview(width = 5, height = 5)
ggsave("results/figures/mouse_cmah.pdf", width = 5, height = 5)
