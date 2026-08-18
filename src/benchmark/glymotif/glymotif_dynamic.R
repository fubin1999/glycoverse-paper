library(glycoverse)
library(tidyverse)
library(glydb)

# Dynamic motifs-----
glycans <- c(
  "Neu5Ac(a2-3)Gal(b1-3)[Neu5Ac(a2-6)]GalNAc(a1-",
  "Gal(b1-4)GlcNAc(b1-6)[Gal(b1-3)]GalNAc(a1-",
  "Gal(b1-3)[Neu5Ac(a2-6)]GalNAc(a1-"
)
motifs <- extract_motif(glycans)

export_cartoons(glycans, "results/figures/glycans_for_dynamic_motifs")
export_cartoons(motifs, "results/figures/dynamic_motifs")

# Dynamic branch motifs-----
glycans <- c(
  "WURCS=2.0/6,20,19/[a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5][a1122h-1a_1-5][a2112h-1b_1-5][Aad21122h-2a_2-6_5*NCC/3=O][a1221m-1a_1-5]/1-1-2-3-1-4-5-1-4-5-3-1-4-5-1-4-1-4-5-6/a4-b1_a6-t1_b4-c1_c3-d1_c6-k1_d2-e1_d4-h1_e4-f1_f3-g2_h4-i1_i3-j2_k2-l1_k6-o1_l4-m1_m3-n2_o4-p1_p3-q1_q4-r1_r3-s2",
  "WURCS=2.0/6,11,10/[a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5][a1122h-1a_1-5][a2112h-1b_1-5][Aad21122h-2a_2-6_5*NCC/3=O][a1221m-1a_1-5]/1-1-2-3-1-4-3-1-4-5-6/a4-b1_a6-k1_b4-c1_c3-d1_c6-g1_d2-e1_e4-f1_g2-h1_h4-i1_i3-j2",
  "WURCS=2.0/6,15,14/[a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5][a1122h-1a_1-5][a2112h-1b_1-5][Aad21122h-2a_2-6_5*NCC/3=O][a1221m-1a_1-5]/1-1-2-3-1-4-5-1-6-4-5-3-1-4-5/a4-b1_b4-c1_c3-d1_c6-l1_d2-e1_d4-h1_e4-f1_f6-g2_h3-i1_h4-j1_j3-k2_l2-m1_m4-n1_n6-o2",
  "WURCS=2.0/8,18,17/[a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5][a1122h-1a_1-5][a2112h-1b_1-5][Aad21122h-2a_2-6_5*NCC/3=O][a2112h-1b_1-5_3*OSO/3=O/3=O][a2112h-1b_1-5_2*NCC/3=O][a1221m-1a_1-5]/1-1-2-3-1-4-5-1-6-3-1-4-5-7-1-4-5-8/a4-b1_a6-r1_b4-c1_c3-d1_c6-j1_d2-e1_d4-h1_e4-f1_f3-g2_h4-i1_j2-k1_j6-o1_k4-l1_l3-m2_l4-n1_o4-p1_p3-q2"
)
motifs <- extract_branch_motif(glycans)
export_cartoons(glycans, "results/figures/glycans_for_dynamic_branch_motifs")
export_cartoons(motifs, "results/figures/dynamic_branch_motifs")

# Dynamic motifs for all O-GalNAc glycans-----
o_glycans <- glydb_structures(glycan_type = "O-GalNAc")  # 731 glycans
count_motif_res <- count_motifs(o_glycans, dynamic_motifs(max_size = 5))  # 1226 motifs
write.csv(count_motif_res, "results/data/glymotif_dynamic_motif_counts.csv")

# Dynamic branch motifs for all N-glycans-----
n_glycans <- glydb_structures(glycan_type = "N")  # 2573 glycans
n_glycans <- n_glycans[have_motif(n_glycans, n_glycan_core(), alignment = "core")]  # 2237 glycans
count_branch_motif_res <- count_motifs(n_glycans, branch_motifs())
write.csv(count_branch_motif_res, "results/data/glymotif_branch_motif_counts.csv")
