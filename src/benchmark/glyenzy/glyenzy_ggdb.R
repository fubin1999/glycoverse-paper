library(tidyverse)
library(glyenzy)
library(glydb)

# prepare data
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

# simulate enzyme reactions
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
  )

write_csv(result_data, "results/data/glyenzy_ggdb_prediction.csv")
