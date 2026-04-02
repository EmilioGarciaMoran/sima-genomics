# FIGURAS PARA JHE SUBMISSION
# Sima de los Huesos — Genomic Identity Study

## 📋 CONTENIDOS

```
figure1_parallel_trees.R       ← Script para árbol mito-nuclear
figure2_heatmap_FINAL.R        ← Script para heatmap distancias nucleares
AJBA_FINAL_DEFENSIBLE.md       ← Manuscrito completo (con captions actualizados)
FIGURES_README.txt             ← Este archivo
```

---

## 🚀 CÓMO GENERAR LAS FIGURAS

### REQUISITOS:

R 4.0+ con paquetes:
```R
install.packages(c("ape", "ggplot2", "reshape2"))
```

### PASOS:

#### 1. FIGURE 1 (Árboles paralelos mito-nuclear)

```R
source("figure1_parallel_trees.R")
```

**Genera:**
- `figure1_parallel_trees.pdf` (vector, para edición)
- `figure1_parallel_trees.tiff` (600 dpi, para submission JHE)
- `figure1_caption.txt` (caption detallado)

**Características:**
- Panel izquierdo: mtDNA rooted (IQ-TREE HKY+G4)
  - Sima basal a (Denisova + Neandertal)
  - Bootstrap 93% annotado
  - Branch lengths a escala (substitutions/site)
  
- Panel derecho: nDNA chromosome 1 (ANGSD distances)
  - Sima-Neandertal CLOSEST (d=0.001376) con red box
  - Branch lengths a escala (pairwise distances)

**Dimensiones:** 14 × 8 inches | 600 dpi

---

#### 2. FIGURE 2 (Heatmap distancias nucleares)

```R
source("figure2_heatmap_FINAL.R")
```

**Genera:**
- `figure2_heatmap_FINAL.pdf` (vector)
- `figure2_heatmap_FINAL.tiff` (600 dpi, para submission)

**Características:**
- Matriz 4×4: Sima, Denisova, Neandertal, Ust'-Ishim
- Dark teal = CLOSEST (0.001376)
- Light teal = FARTHEST (0.004643)
- Red border: Sima-Neandertal minimum distance
- Valores numéricos anotados en cada celda

**Dimensiones:** 7.5 × 6.0 inches | 600 dpi

---

## 📝 CAPTIONS (integrados en manuscrito)

### FIGURE 1 CAPTION:

"**Figure 1. Parallel Mitochondrial and Nuclear Phylogenies Reveal Mito-Nuclear Discordance 
in Sima de los Huesos.** (Left panel) Rooted mitochondrial phylogeny (IQ-TREE v3.0.1, HKY+G4 
model, 1,000 ultrafast bootstrap replicates). Sima de los Huesos mtDNA occupies a basal sister 
position to the clade containing Denisova and Neandertal specimens (bootstrap support 93%). 
Tree explicitly rooted on Pan troglodytes (outgroup, branch length truncated). Branch lengths 
proportional to substitutions per site. Maximum likelihood score: −30678.608; Robinson-Foulds 
distance from consensus tree = 0. (Right panel) Nuclear phylogeny based on pairwise genetic 
distances (ANGSD v0.941, chromosome 1, ~2.19 million filtered sites). The minimum pairwise 
distance is Sima–Neandertal (d = 0.001376, red border), approximately 1.82-fold closer than 
Sima–Denisova (d = 0.002506). Outgroup: Ust'-Ishim (45 ka, modern human). Branch lengths 
proportional to proportion of discordant genotype calls. The apparent discordance—Sima mtDNA 
basal but nuclear genome closest to Neandertal—is consistent with incomplete lineage sorting 
in a transitional population with deep Middle Pleistocene roots. This pattern reflects distinct 
coalescence times for the maternal and autosomal genomes."

### FIGURE 2 CAPTION:

"**Figure 2. Pairwise Nuclear Genetic Distance Matrix (Chromosome 1).** Colour-coded heatmap 
of pairwise genetic distances computed from ~2.19 million chromosome 1 sites (ANGSD v0.941 
with post-mortem damage filtering and block jackknife error correction). Darker teal indicates 
closer genetic relatedness; lighter shades indicate greater divergence. Critical feature: The 
minimum pairwise distance is Sima–Neandertal (d = 0.001376, indicated by red border), 
approximately 1.82-fold shorter than Sima–Denisova distance (d = 0.002506). This pattern 
demonstrates unambiguous preferential nuclear affinity between Sima and Neandertals, contrasting 
with mtDNA topology where Sima is basal to both lineages. Diagonal entries marked "—" (same 
taxon). Reference outgroup: Ust'-Ishim (45 ka, modern human). Data: ANGSD genotype likelihoods 
aligned to hg19."

---

## ✅ JHE SUBMISSION CHECKLIST

- [x] Figure 1: Two parallel trees (mito | nDNA) with branch lengths
- [x] Figure 2: Heatmap with Sima-Neandertal highlighted (minimum distance)
- [x] Captions: Detailed, explaining mito-nuclear discordance
- [x] Resolution: 600 dpi TIFF for both figures
- [x] Format: PDF (vector) + TIFF (raster)
- [x] Bootstrap/Z-scores: In captions (not on figures)

---

## 📊 DATOS CLAVE (para verificar figuras)

### mtDNA Tree (IQ-TREE HKY+G4):
- LogL: −30678.608
- Bootstrap: 93% (Sima basal)
- Topology: (Pan, (Sima, (Denisova, (Neandertal, Feldhofer, Mezmaiskaya))))

### nDNA Distances (ANGSD chr1):
- Sima-Neandertal: 0.001376 ← MINIMUM
- Sima-Denisova: 0.002506
- Ratio: 1.82× (Sima closer to Neandertal)

---

## 🔍 TROUBLESHOOTING

**Error: "No such file or directory"**
→ Asegúrate de estar en el directorio correcto donde están los scripts R

**Error: "package 'ape' is required"**
→ Instala: `install.packages("ape")`

**Figuras se ven en blanco/distorsionadas**
→ Verifica resolución: TIFF debe ser 600 dpi mínimo

**Problema con fuentes (Helvetica/sans)**
→ Usa `extrafont::loadfonts()` o cambia a `family = "sans"` (por defecto)

---

## 📮 PARA ENVÍO A JHE

1. Usa los archivos `.tiff` (600 dpi) para submission
2. Guarda los `.pdf` como backup (vector, editable)
3. Captions van en el archivo Word/manuscript, NO en las figuras
4. Formato JHE: Figuras separadas (no incrustadas en texto)

---

## HISTORIAL DE CAMBIOS

**v2.0 (Final):**
- Figure 1: Dos árboles paralelos (reemplaza tanglegram)
- Figure 2: Heatmap optimizado con red border en mínimo
- Captions reescritos para claridad mito-nuclear discordance
- 600 dpi TIFF para JHE requirement
