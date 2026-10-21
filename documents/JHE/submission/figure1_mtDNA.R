# OBSOLETO (2026-10-21): NO USAR. Ver CORRECTION.md, sección "Figuras".
# El árbol está escrito a mano (Newick literal): usa longitudes de rama reales de IQ-TREE sobre una topología 'Sima basal' que no es la que dio el análisis, y el '93%' se coloca con text() en coordenadas fijas.
# Sustituido por mtDNA_repro/scripts/06_figure.py, que dibuja los árboles leyendo la salida de IQ-TREE.
# ============================================================
# Figure 1 — Mitochondrial DNA phylogeny ONLY
# IQ-TREE v3.0.1 · HKY+G4 · rooted on Pan troglodytes
# Output: figure1_mtDNA.pdf + figure1_mtDNA.tiff (600 dpi)
# ============================================================

library(ape)

mt_tree <- read.tree(text =
  "(Pan_troglodytes:0.1320,(Sima_430ka:0.0048,(Denisova3:0.0097,(Neandertal_Vindija:0.0003,Feldhofer:0.0003):0.0014):0.0089):0.0052);")

mt_tree$tip.label <- c(
  "Pan troglodytes",
  "Sima de los Huesos  (430 ka)",
  "Denisova  (70 ka)",
  "H. neanderthalensis  Vindija  (44 ka)",
  "H. neanderthalensis  Feldhofer  (40 ka)"
)

# Edge colours: Sima branch in red
edge_col <- rep("black", nrow(mt_tree$edge))
sima_tip  <- which(grepl("Sima", mt_tree$tip.label))
for (i in seq_len(nrow(mt_tree$edge)))
  if (mt_tree$edge[i, 2] == sima_tip) edge_col[i] <- "#C94C4C"

draw_fig1 <- function() {
  par(mar = c(3, 1, 2, 8), mgp = c(2, 0.6, 0),
      family = "sans")

  plot(mt_tree,
       type            = "phylogram",
       use.edge.length = TRUE,
       show.tip.label  = TRUE,
       edge.color      = edge_col,
       edge.width      = 2.8,
       cex             = 1.20,
       label.offset    = 0.002,
       x.lim           = c(0, 0.18))

  # Scale bar
  sw <- 0.03; sx <- 0.003; sy <- 0.55
  rect(sx, sy - 0.07, sx + sw, sy + 0.07,
       col = "black", border = NA)
  text((2*sx + sw)/2, sy - 0.32,
       "0.03 substitutions/site", cex = 1.0, adj = 0.5)

  # Bootstrap value on key node
  text(0.005, 4.60, "93%", cex = 1.6, font = 2, col = "#C94C4C")

  # Outgroup label
  mtext(expression(italic("Pan troglodytes") ~ "(outgroup)"),
        side = 1, line = 0.5, adj = 0, cex = 0.80, col = "grey45")
}

# PDF
pdf("figure1_mtDNA.pdf", width = 10, height = 7)
draw_fig1()
dev.off()

# TIFF 600 dpi
tiff("figure1_mtDNA.tiff",
     width = 10, height = 7, units = "in",
     res = 600, compression = "lzw")
draw_fig1()
dev.off()

message("figure1_mtDNA.pdf + .tiff saved")
