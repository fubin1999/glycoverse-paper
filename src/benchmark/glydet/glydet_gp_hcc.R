library(glycoverse)
library(tidyverse)
library(ggrepel)
library(patchwork)

exp <- auto_clean(real_experiment) |>
  filter_col(group %in% c("H", "Y"), .drop_levels = TRUE)

motifs <- c(
  lewis_by = "dHex(??-?)Hex(??-?)[dHex(??-?)]HexNAc(??-",
  lewis_ax = "Hex(??-?)[dHex(??-?)]HexNAc(??-",
  sia_lewis_ax = "NeuAc(??-?)Hex(??-?)[dHex(??-?)]HexNAc(??-"
)
lewis_exp <- quantify_motifs(exp, motifs, method = "absolute")
lewis_dea_res <- lewis_exp |>
  gly_ttest() |>
  get_tidy_result()

lewis_dea_res |> filter(p_adj < 0.05)

plot_volcano <- function(dea_res) {
  dea_res |>
    mutate(
      regulate = case_when(
        p_adj < 0.05 & log2fc > 1 ~ "up",
        p_adj < 0.05 & log2fc < 1 ~ "down",
        .default = "none"
      ),
      site = if_else(
        regulate != "none",
        paste0(gene, "-N", protein_site),
        NA_character_
      )
    ) |>
    ggplot(aes(log2fc, -log10(p_adj))) +
    geom_point(aes(color = regulate, fill = regulate), size = 4) +
    geom_text_repel(aes(label = site, color = regulate)) +
    scale_color_manual(values = c(up = "#D2755B", down = "#818CA5", none = "grey90")) +
    scale_fill_manual(values = c(up = "#D2755B", down = "#818CA5", none = "grey90")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    guides(color = "none", fill = "none") +
    theme_bw()
}

volcano_lewis_ax <- lewis_dea_res |>
  filter(trait == "lewis_ax", !is.na(p_val)) |>
  plot_volcano()

volcano_lewis_by <- lewis_dea_res |>
  filter(trait == "lewis_by", !is.na(p_val)) |>
  plot_volcano()

volcano_sia_lewis_ax <- lewis_dea_res |>
  filter(trait == "sia_lewis_ax", !is.na(p_val)) |>
  plot_volcano()

volcanos <- volcano_lewis_ax + volcano_lewis_by + volcano_sia_lewis_ax
ggsave("results/figures/lewis_volcanos.pdf", plot = volcanos, width = 12, height = 4)
