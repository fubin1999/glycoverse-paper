library(tidyverse)
library(glycoverse)

fully_determined <- read_csv("data/benchmark/text_nomenclatures/glycan_fully_determined.csv")
glycoct <- read_csv("data/benchmark/text_nomenclatures/glycan_sequences_glycoct.csv") |>
  mutate(sequence_glycoct = str_replace_all(sequence_glycoct, " ", "\n"))
wurcs <- read_csv("data/benchmark/text_nomenclatures/glycan_sequences_wurcs.csv")
iupac_ext <- read_csv("data/benchmark/text_nomenclatures/glycan_sequences_iupac_extended.csv")
iupac_con <- read_csv("data/benchmark/text_nomenclatures/glycan_sequences_iupac_condensed.csv")

glycoct2 <- fully_determined |>
  left_join(glycoct) |>
  filter(!str_detect(sequence_glycoct, "aldi"))
glycoct2_glycans <- parse_glycoct(glycoct2$sequence_glycoct, on_failure = "na")
sum(!is.na(glycoct2_glycans)) / length(glycoct2_glycans)

wurcs2 <- wurcs |>
  right_join(glycoct2)
wurcs2_glycans <- parse_wurcs(wurcs2$sequence_wurcs, on_failure = "na")
sum(!is.na(wurcs2_glycans)) / length(wurcs2_glycans)
                                     