library(tidyverse)
library(glycoverse)
library(bench)

accession <- read_csv("data/text_nomenclatures/glycan_fully_determined.csv")
glycoct <- read_csv("data/text_nomenclatures/glycan_sequences_glycoct.csv")
iupac <- read_csv("data/text_nomenclatures/glycan_sequences_iupac_extended.csv")
wurcs <- read_csv("data/text_nomenclatures/glycan_sequences_wurcs.csv")

glycoct <- glycoct |>
  semi_join(accession, by = "glytoucan_ac") |>
  rename(sequence = sequence_glycoct)

iupac <- iupac |>
  semi_join(accession, by = "glytoucan_ac") |>
  rename(sequence = sequence_iupac_extended)

wurcs <- wurcs |>
  semi_join(accession, by = "glytoucan_ac") |>
  rename(sequence = sequence_wurcs)

run_parsing <- function(type, n) {
  checkmate::assert_choice(type, c("glycoct", "iupac", "wurcs"))
  checkmate::assert_int(n)
  data <- switch(type, glycoct = glycoct, iupac = iupac, wurcs = wurcs)
  strings <- sample(data$sequence, n, replace = FALSE)
  tryCatch(auto_parse(strings), error = function(e) NULL)
}

results <- bench::press(
  type = c("glycoct", "iupac", "wurcs"),
  n = c(100, 200, 300, 400, 500),
  bench::mark(run_parsing(type, n))
)

ggplot(results, aes(n, median, color = type)) +
  geom_line() +
  geom_point() +
  scale_color_viridis_d(labels = c("GlycoCT", "IUPAC", "WURCS")) +
  labs(x = "Number of Sequences", y = "Median Time (s)", color = "Type") +
  theme_minimal()
ggpreview(width = 6, height = 4)
ggsave("results/figures/glyparse_bench.pdf", width = 5, height = 4)
