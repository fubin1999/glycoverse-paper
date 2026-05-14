library(glycoverse)
library(tidyverse)

exp1 <- read_byonic_pglycoquant("data/glycoproteomics/byonic_pglycoquant.list")
df1 <- as_tibble(exp)
write_csv(df1, "results/data/byonic_pglycoquant.csv")

exp2 <- read_byonic_byologic("data/glycoproteomics/byonic_byologic.csv")
df2 <- as_tibble(exp2)
write_csv(df2, "results/data/byonic_byologic.csv")
