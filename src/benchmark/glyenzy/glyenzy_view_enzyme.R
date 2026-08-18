library(glyrepr)
library(glyenzy)
library(purrr)

glycan <- "Neu5Ac(a2-3)Gal(b1-3)[Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-6)]GalNAc(a1-"
enzymes <- c(
  "GALNT1", "C1GALT1", "GCNT1", "ST3GAL1", "ST3GAL3",
  "B4GALT1", "FUT3", "FUT6"
)

cartoons <- map(enzymes, ~ view_enzyme(glycan, .x))
dirpath <- "results/figures/view_enzyme_cartoons"
fs::dir_create(dirpath)
filenames <- map_chr(enzymes, ~ fs::path(dirpath, paste0(.x, ".pdf")))
walk2(cartoons, filenames, save_cartoon)
