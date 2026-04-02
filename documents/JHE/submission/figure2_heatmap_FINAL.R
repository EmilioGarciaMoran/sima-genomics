# Figure 2: Pairwise Nuclear Genetic Distance Heatmap
# Sima de los Huesos — Journal of Human Evolution
# Data: chromosome 1, ~2.19M subsampled sites, ANGSD genotype likelihoods
# R packages required: ggplot2, reshape2

# install.packages(c("ggplot2", "reshape2"))
library(ggplot2)
library(reshape2)

# ────────────────────────────────────────────────────────────────────────────
# DISTANCE MATRIX (proportion of discordant genotype calls)
# ────────────────────────────────────────────────────────────────────────────

taxa <- c("Sima\n~430 ka", "Denisova\n~70 ka", "Neandertal\n~50 ka", "Ust'-Ishim\n~45 ka")

mat <- matrix(
  c(NA,       0.002506, 0.001376, 0.002557,
    0.002506, NA,       0.003481, 0.004643,
    0.001376, 0.003481, NA,       0.003496,
    0.002557, 0.004643, 0.003496, NA),
  nrow = 4, byrow = TRUE,
  dimnames = list(taxa, taxa)
)

# ────────────────────────────────────────────────────────────────────────────
# MELT TO LONG FORMAT
# ────────────────────────────────────────────────────────────────────────────

df <- melt(mat, varnames = c("Taxon1", "Taxon2"), value.name = "Distance")
df$Taxon1 <- factor(df$Taxon1, levels = taxa)
df$Taxon2 <- factor(df$Taxon2, levels = rev(taxa))

# Label: show value or "—" for diagonal
df$label <- ifelse(is.na(df$Distance), "—", sprintf("%.4f", df$Distance))

# Highlight minimum off-diagonal (Sima–Neandertal = CLOSEST)
df$is_min <- !is.na(df$Distance) & df$Distance == min(df$Distance, na.rm = TRUE)

# ────────────────────────────────────────────────────────────────────────────
# COLOR PALETTE
# Dark teal = CLOSEST (0.001376)
# Light = FARTHEST (0.004643)
# ────────────────────────────────────────────────────────────────────────────

pal_low  <- "#0d3b47"   # very dark teal (closest)
pal_mid  <- "#5b9ab5"   # mid blue
pal_high <- "#f0f4f8"   # near-white (most divergent)

# ────────────────────────────────────────────────────────────────────────────
# PLOT WITH GGPLOT2
# ────────────────────────────────────────────────────────────────────────────

p <- ggplot(df, aes(x = Taxon1, y = Taxon2, fill = Distance)) +
  
  # Main tiles
  geom_tile(color = "white", linewidth = 1.0) +
  
  # Text labels
  geom_text(aes(label = label,
                color = ifelse(is.na(Distance) | Distance > 0.003, "dark", "light")),
            size = 4.2, fontface = "plain", family = "sans") +
  
  # RED BORDER on minimum cell (Sima-Neandertal)
  geom_tile(data = subset(df, is_min),
            aes(x = Taxon1, y = Taxon2),
            fill = NA, color = "#c0392b", linewidth = 2.2) +
  
  # Color gradient
  scale_fill_gradient2(
    low      = pal_low,
    mid      = pal_mid,
    high     = pal_high,
    midpoint = 0.003,
    na.value = "#e8e8e8",
    name     = "Pairwise\ndistance",
    limits   = c(0.001, 0.005),
    breaks   = c(0.001, 0.002, 0.003, 0.004, 0.005),
    labels   = c("0.001\n(closest)", "0.002", "0.003", "0.004", "0.005\n(most divergent)")
  ) +
  
  # Text color
  scale_color_manual(values = c("dark" = "grey20", "light" = "white"),
                     guide  = "none") +
  
  # Layout
  scale_x_discrete(position = "top") +
  
  # Labels
  labs(
    title    = "Pairwise Nuclear Genetic Distances",
    subtitle = "Chromosome 1 · ~2.19 million sites · proportion of discordant genotype calls",
    caption  = "Red border: minimum pairwise distance (Sima–Neandertal, d = 0.001376).\nDiagonal not applicable (same taxon). Data: ANGSD genotype likelihoods, hg19."
  ) +
  
  # Theme
  theme_minimal(base_family = "sans", base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 14, hjust = 0,
                                    margin = margin(b = 5)),
    plot.subtitle    = element_text(size = 10, color = "grey40", hjust = 0,
                                    margin = margin(b = 15)),
    plot.caption     = element_text(size = 9, color = "grey50", hjust = 0,
                                    margin = margin(t = 12), lineheight = 1.3),
    axis.text.x      = element_text(size = 11, face = "bold", color = "grey20",
                                    hjust = 0.5),
    axis.text.y      = element_text(size = 11, face = "bold", color = "grey20",
                                    hjust = 1),
    axis.title       = element_blank(),
    panel.grid       = element_blank(),
    legend.title     = element_text(size = 10, color = "grey30", face = "plain"),
    legend.text      = element_text(size = 9, color = "grey30"),
    legend.key.height = unit(1.4, "cm"),
    plot.margin      = margin(15, 15, 15, 15),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

# ────────────────────────────────────────────────────────────────────────────
# EXPORT: PDF (vector)
# ────────────────────────────────────────────────────────────────────────────

ggsave("figure2_heatmap_FINAL.pdf",
       plot   = p,
       width  = 7.5,
       height = 6.0,
       units  = "in",
       device = cairo_pdf,
       dpi    = 300)

# ────────────────────────────────────────────────────────────────────────────
# EXPORT: TIFF 600 dpi (JHE requirement)
# ────────────────────────────────────────────────────────────────────────────

ggsave("figure2_heatmap_FINAL.tiff",
       plot   = p,
       width  = 7.5,
       height = 6.0,
       units  = "in",
       dpi    = 600,
       compression = "lzw")

message("✓ figure2_heatmap_FINAL.pdf saved (vector)")
message("✓ figure2_heatmap_FINAL.tiff saved (600 dpi, ready for JHE submission)")
