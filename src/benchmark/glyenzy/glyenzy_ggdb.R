library(tidyverse)
library(glyenzy)
library(glydb)

# prepare data-----
data <- read_csv("data/benchmark/glycosmos_ggdbs_pubmed.csv") |>
  janitor::clean_names() |>
  filter(reaction != "General_reaction") |>
  mutate(
    acceptor_structure = glytoucan_to_struc(acceptor),
    product_structure = glytoucan_to_struc(product)
  ) |>
  filter(!is.na(acceptor_structure)) |>
  filter(!(reaction == "Substrate_transferred" & is.na(product_structure))) |>
  mutate(gene_symbol = if_else(str_detect(gene_symbol, "ABO"), "ABO", gene_symbol)) |>
  filter(gene_symbol %in% db_enzymes(return_str = TRUE))

# simulate enzyme reactions-----
result_data <- data |>
  rowwise() |>
  mutate(
    predict_products = list(apply_enzyme(acceptor_structure, gene_symbol)),
    matched = if_else(
      reaction == "Substrate_transferred",
      product_structure %in% predict_products,
      length(predict_products) == 0
    ),
    predict_products = paste(predict_products, collapse = ";")
  ) |>
  unnest()

write_csv(result_data, "results/data/glyenzy_ggdb_prediction.csv")

# plot-----
p <- result_data |>
  summarise(n = n(), .by = c(reaction, matched)) |>
  ggplot(aes(n, reaction)) +
  geom_col(aes(fill = matched), position = position_stack(), width = 0.6) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c("#f9c22b", "#79add2")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  annotate("text", x = 5, y = 2.45, label = "Substrates transfered", hjust = 0) +
  annotate("text", x = 5, y = 1.45, label = "Substrates not transfered", hjust = 0) +
  theme_classic() +
  theme(
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank()
  )

ggimage::ggpreview(plot = p, width = 4, height = 1.5)
ggsave("results/figures/glyenzy_ggdb.pdf", plot = p, width = 4, height = 1.5)
