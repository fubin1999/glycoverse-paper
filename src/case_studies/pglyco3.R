library(tidyverse)
library(glycoverse)

# Data reading and preprocessing----
exp <- read_pglyco3("data/glycoproteomics/pglyco3_result.txt", parse_structure = TRUE)
sample_info <- read_csv("data/glycoproteomics/sample_info.csv") |>
  mutate(
    sample = str_split_i(Name, fixed("."), 1),
    protease = case_when(
      Trypsin == "+" ~ "Trypsin",
      `Trypsin+Glu-C` == "+" ~ "Trypsin+Glu-C",
      `Trypsin+Lysine-C` == "+" ~ "Trypsin+Lysine-C",
      .default = "others"
    ),
    lc = case_when(
      `Sepharose CL-4B` == "+" ~ "Sepharose CL-4B",
      `ZIC-HILIC` == "+" ~ "ZIC-HILIC",
      .default = "others"
    ),
    organ = case_when(
      Heart == "+" ~ "heart",
      liver == "+" ~ "liver",
      Brain == "+" ~ "brain",
      lung == "+" ~ "lung",
      kidney == "+" ~ "kidney",
      .default = "others"
    ),
    condition = str_c(organ, "_", protease, "_", lc),
    protease = factor(protease),
    lc = factor(lc),
    organ = factor(organ),
    condition = factor(condition),
    .keep = "unused"
  )

clean_exp <- exp |>
  left_join_obs(sample_info, by = "sample") |>
  normalize_median() |>
  impute_zero() |>
  aggregate(to_level = "gfs")

write_rds(clean_exp, "results/data/glycoproteomics_clean_exp.rds")
clean_exp <- read_rds("results/data/glycoproteomics_clean_exp.rds")

# Derived trait analysis-----
trait_exp <- derive_traits(clean_exp)
trait_anova_res <- gly_anova(trait_exp, group_col = "organ")
get_tidy_result(trait_anova_res) |>
  filter(p_adj < 0.05)