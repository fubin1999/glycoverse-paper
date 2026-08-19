library(glyenzy)
library(glyrepr)
library(glydraw)
library(igraph)

glycan <- "Neu5Ac(a2-3)Gal(b1-3)[Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-6)]GalNAc(a1-"
path <- trace_biosynthesis(glycan)
class(path) <- "igraph"

pdf("results/figures/glyenzy_path.pdf", width = 6, height = 6)
plot(
  path,
  layout = layout_as_tree(path),
  vertex.size = 4,
  vertex.frame.color = "#abc6cf",
  vertex.frame.width = 3,
  vertex.color = "white",
  edge.arrow.size = 0.3,
  margin = 0
)
dev.off()

middle_glycans <- V(path)$name
export_cartoons(
  middle_glycans,
  "results/figures/glyenzy_path_middle_glycans/",
  file_ext = "pdf"
)

glycan2 <- "Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc(b1-"
path2 <- trace_biosynthesis(glycan2)

pdf("results/figures/glyenzy_path2.pdf", width = 6, height = 6)
plot(
  path2,
  layout = layout_as_tree(path2),
  vertex.size = 4,
  vertex.frame.color = "#abc6cf",
  vertex.frame.width = 3,
  vertex.color = "white",
  edge.arrow.size = 0.3,
  margin = 0
)
dev.off()

middle_glycans2 <- V(path2)$name
export_cartoons(
  middle_glycans2,
  "results/figures/glyenzy_path2_middle_glycans/",
  file_ext = "pdf"
)
