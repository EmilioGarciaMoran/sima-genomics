# Figure 2: Pairwise nuclear genetic distance heatmap
# Sima de los Huesos — Journal of Human Evolution
# Data: chromosome 1, ~2.19M subsampled sites, ANGSD genotype likelihoods
# R packages required: ggplot2, reshape2, ggtext

# install.packages(c("ggplot2", "reshape2", "ggtext"))
library(ggplot2)
library(reshape2)

# ── Distance matrix (proportion of discordant genotype calls) ──────────────
taxa <- c("Sima\n~430 ka", "Denisova\n~70 ka", "Neandertal\n~50 ka", "Ust'-Ishim\n~45 ka")

mat <- matrix(
  c(NA,       0.002506, 0.001376, 0.002557,
    0.002506, NA,       0.003481, 0.004643,
    0.001376, 0.003481, NA,       0.003496,
    0.002557, 0.004643, 0.003496, NA),
  nrow = 4, byrow = TRUE,
  dimnames = list(taxa, taxa)
)

# ── Melt to long format ────────────────────────────────────────────────────
df <- melt(mat, varnames = c("Taxon1", "Taxon2"), value.name = "Distance")
df$Taxon1 <- factor(df$Taxon1, levels = taxa)
df$Taxon2 <- factor(df$Taxon2, levels = rev(taxa))  # reverse y for readability

# ── Label: show value or "—" for diagonal ─────────────────────────────────
df$label <- ifelse(is.na(df$Distance), "—", sprintf("%.4f", df$Distance))

# ── Highlight the minimum off-diagonal (Sima–Neandertal) ──────────────────
df$is_min <- !is.na(df$Distance) & df$Distance == min(df$Distance, na.rm = TRUE)

# ── Color palette: white → steel blue (low = close = warm, high = far = cool)
# Reversed: closer = darker teal, farther = lighter
pal_low  <- "#1a3a4a"   # dark teal  → closest
pal_mid  <- "#5b9ab5"   # mid blue
pal_high <- "#f0f4f8"   # near-white → most divergent

# ── Plot ───────────────────────────────────────────────────────────────────
p <- ggplot(df, aes(x = Taxon1, y = Taxon2, fill = Distance)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = label,
                color = ifelse(is.na(Distance) | Distance > 0.003, "dark", "light")),
            size = 3.8, fontface = "plain", family = "Helvetica") +
  # Red border on minimum cell
  geom_tile(data = subset(df, is_min),
            aes(x = Taxon1, y = Taxon2),
            fill = NA, color = "#c0392b", linewidth = 1.8) +
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
  scale_color_manual(values = c("dark" = "grey20", "light" = "white"),
                     guide  = "none") +
  scale_x_discrete(position = "top") +
  labs(
    title    = "Pairwise Nuclear Genetic Distances",
    subtitle = "Chromosome 1 · ~2.19 million sites · proportion of discordant genotype calls",
    caption  = "Red border: minimum pairwise distance (Sima–Neandertal, d = 0.001376).\nDiagonal not applicable (same taxon). Data: ANGSD genotype likelihoods, hg19."
  ) +
  theme_minimal(base_family = "Helvetica", base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 13, hjust = 0,
                                    margin = margin(b = 4)),
    plot.subtitle    = element_text(size = 9, color = "grey40", hjust = 0,
                                    margin = margin(b = 12)),
    plot.caption     = element_text(size = 8, color = "grey50", hjust = 0,
                                    margin = margin(t = 10)),
    axis.text.x      = element_text(size = 10, face = "bold", color = "grey20",
                                    hjust = 0.5),
    axis.text.y      = element_text(size = 10, face = "bold", color = "grey20",
                                    hjust = 1),
    axis.title       = element_blank(),
    panel.grid       = element_blank(),
    legend.title     = element_text(size = 9, color = "grey30"),
    legend.text      = element_text(size = 8, color = "grey30"),
    legend.key.height = unit(1.2, "cm"),
    plot.margin      = margin(20, 20, 20, 20),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

# ── Export ─────────────────────────────────────────────────────────────────
# PDF (vector, for journal submission)
ggsave("figure2_heatmap.pdf",
       plot   = p,
       width  = 7,
       height = 5.5,
       units  = "in",
       device = cairo_pdf)

# TIFF 600 dpi (JHE requires minimum 300 dpi for halftones, 600 for figures)
ggsave("figure2_heatmap.tiff",
       plot   = p,
       width  = 7,
       height = 5.5,
       units  = "in",
       dpi    = 600,
       compression = "lzw")

message("✓ figure2_heatmap.pdf and .tiff saved")
