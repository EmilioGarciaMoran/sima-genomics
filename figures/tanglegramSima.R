# Figure 1: Mito-nuclear tanglegram — Sima de los Huesos
# Journal of Human Evolution submission
# R packages required: ape, dendextend

# install.packages(c("ape", "dendextend"))
library(ape)
library(dendextend)

# ── Trees — ingroup only (outgroups omitted, noted in caption) ─────────────
# Mitochondrial (Phase 2): Sima basal sister to (Denisova + Neandertal)
mt_tree <- read.tree(text =
  "(Sima_430ka:7,(Denisova_70ka:3,Neandertal_50ka:3):4);")

# Nuclear (Phase 3): Sima sister to Neandertal; Denisova outgroup of pair
nuc_tree <- read.tree(text =
  "(Denisova_70ka:5,(Sima_430ka:3,Neandertal_50ka:3):2);")

# ── Convert to dendrogram via cophenetic distances ─────────────────────────
d1 <- as.dendrogram(hclust(as.dist(cophenetic(mt_tree)), method = "average"))
d2 <- as.dendrogram(hclust(as.dist(cophenetic(nuc_tree)), method = "average"))

# ── Export ─────────────────────────────────────────────────────────────────
pdf("figure1_tanglegram.pdf", width = 11, height = 7)

par(oma = c(0, 0, 2.5, 0))

tanglegram(
  d1, d2,
  highlight_distinct_edges    = FALSE,
  common_subtrees_color_lines = FALSE,
  margin_inner   = 14,
  margin_outer   = 3,
  lab.cex        = 0.85,
  lwd            = 2,
  edge.lwd       = 2,
  columns_width  = c(5, 2, 5),
  color_lines    = c("#C9A84C", "#C94C4C", "#4CB8C9", "#8B7A6B"),
  main_left      = "Mitochondrial DNA (Phase 2)",
  main_right     = "Nuclear DNA — chr1 (Phase 3)",
  main           = "",
  cex_main       = 0.85
)
mtext("Sima de los Huesos — Mito-Nuclear Discordance",
      side = 3, outer = TRUE, cex = 1.1, font = 2, line = 0.8)

dev.off()

# ── TIFF 600 dpi for JHE submission ───────────────────────────────────────
tiff("figure1_tanglegram.tiff",
     width = 11, height = 7, units = "in",
     res = 600, compression = "lzw")

par(oma = c(0, 0, 2.5, 0))

tanglegram(
  d1, d2,
  highlight_distinct_edges    = FALSE,
  common_subtrees_color_lines = FALSE,
  margin_inner   = 14,
  margin_outer   = 3,
  lab.cex        = 0.85,
  lwd            = 2,
  edge.lwd       = 2,
  columns_width  = c(5, 2, 5),
  color_lines    = c("#C9A84C", "#C94C4C", "#4CB8C9", "#8B7A6B"),
  main_left      = "Mitochondrial DNA (Phase 2)",
  main_right     = "Nuclear DNA — chr1 (Phase 3)",
  main           = "",
  cex_main       = 0.85
)
mtext("Sima de los Huesos — Mito-Nuclear Discordance",
      side = 3, outer = TRUE, cex = 1.1, font = 2, line = 0.8)

dev.off()

# ── Note on Outgroup labels ────────────────────────────────────────────────
# "Outgroup" = Pan troglodytes in the mtDNA tree (left panel)
# "Outgroup" = Ust'-Ishim in the nuclear tree (right panel)
# This asymmetry should be noted in the Figure 1 caption (already in manuscript)

message("✓ figure1_tanglegram.pdf and .tiff saved")
