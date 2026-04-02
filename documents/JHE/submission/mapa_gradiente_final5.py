#!/usr/bin/env python3
# =============================================================
# MAPA GRADIENTE GENÉTICO - ADN NUCLEAR  · FINAL v5
# Figura limpia para publicación en inglés
# Sin texto interpretativo — va en el caption del artículo
# =============================================================

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import geopandas as gpd
from shapely import prepare, contains_xy
import warnings
warnings.filterwarnings("ignore")

# --- 1. SITIOS ----------------------------------------------
sitios = {
    "Sima de los Huesos\n430 ka": {
        "lon": -3.5194, "lat": 42.35,
        "color": np.array([205,  18,  42]) / 255.0
    },
    "Vindija\n44 ka": {
        "lon": 16.084,  "lat": 46.301,
        "color": np.array([228,  62,  12]) / 255.0
    },
    "Ust'-Ishim\n45 ka": {
        "lon": 71.1667, "lat": 57.6833,
        "color": np.array([ 18,  55, 210]) / 255.0
    },
    "Denisova\n70 ka": {
        "lon": 84.6767, "lat": 51.3975,
        "color": np.array([ 55,  35, 180]) / 255.0
    },
}

# --- 2. MAPA BASE -------------------------------------------
print("Loading map...")
world   = gpd.read_file(
    "https://naturalearth.s3.amazonaws.com/110m_cultural/"
    "ne_110m_admin_0_countries.zip"
)
eurasia = world.cx[-15:145, 27:73]
land    = eurasia.unary_union
prepare(land)

coast = gpd.read_file(
    "https://naturalearth.s3.amazonaws.com/110m_physical/"
    "ne_110m_coastline.zip"
).cx[-15:145, 27:73]

# --- 3. GRID IDW --------------------------------------------
res  = 0.20
lons = np.arange(-15,  145 + res, res)
lats = np.arange( 27,   73 + res, res)
LON, LAT = np.meshgrid(lons, lats)
pts_lon = LON.ravel()
pts_lat = LAT.ravel()

print("Filtering land points...")
mask = contains_xy(land, pts_lon, pts_lat)
ll, lt = pts_lon[mask], pts_lat[mask]
print(f"  {mask.sum()} land points")

def idw(ll, lt, channel, pow_=2.0):
    out, w_sum = np.zeros(len(ll)), np.zeros(len(ll))
    for s in sitios.values():
        d = np.hypot(ll - s["lon"], lt - s["lat"])
        w = 1.0 / (d**pow_ + 0.5)
        out   += w * s["color"][channel]
        w_sum += w
    return np.clip(out / w_sum, 0, 1)

print("Computing IDW...")
R = idw(ll, lt, 0)
G = idw(ll, lt, 1)
B = idw(ll, lt, 2)

nrows, ncols = LON.shape
RGBA = np.full((nrows, ncols, 4), np.nan)
ri, ci = np.where(mask.reshape(nrows, ncols))
for k, (i, j) in enumerate(zip(ri, ci)):
    RGBA[i, j] = [R[k], G[k], B[k], 0.68]

# --- 4. FIGURA ----------------------------------------------
print("Rendering...")
fig, ax = plt.subplots(figsize=(180/25.4, 110/25.4),
                       dpi=600, facecolor="white")
ax.set_facecolor("white")

eurasia.plot(ax=ax, color="#F3F1ED", edgecolor="none", zorder=1)
ax.imshow(RGBA, origin="lower", extent=[-15, 145, 27, 73],
          aspect="auto", interpolation="bilinear", zorder=2)
eurasia.plot(ax=ax, color="none", edgecolor="#AAAAAA",
             linewidth=0.18, zorder=3)
coast.plot(ax=ax, color="#2a2a2a", linewidth=0.50, zorder=4)

# --- 5. PUNTOS Y ETIQUETAS (en inglés) ----------------------
offsets = {
    "Sima de los Huesos\n430 ka": ( 6.0, -5.0),
    "Vindija\n44 ka":             ( 3.0,  3.5),
    "Ust'-Ishim\n45 ka":          ( 3.5,  3.8),
    "Denisova\n70 ka":            ( 4.5, -5.5),
}
for label, s in sitios.items():
    ax.plot(s["lon"], s["lat"], "o",
            markersize=5.5, markerfacecolor="white",
            markeredgecolor="#111111", markeredgewidth=0.75, zorder=6)
    dx, dy = offsets[label]
    ax.annotate(label,
        xy=(s["lon"], s["lat"]),
        xytext=(s["lon"]+dx, s["lat"]+dy),
        fontsize=6.5, ha="center", va="center",
        color="#111111", linespacing=1.3,
        bbox=dict(boxstyle="round,pad=0.28", facecolor="white",
                  edgecolor="#BBBBBB", alpha=0.92, linewidth=0.45),
        arrowprops=dict(arrowstyle="-", color="#777777", linewidth=0.45),
        zorder=7)

# --- 6. EJES ------------------------------------------------
ax.set_xlim(-15, 145); ax.set_ylim(27, 73)
ax.set_xticks(range(-10, 141, 20))
ax.set_xticklabels(
    [f"{abs(x)}°{'W' if x<0 else 'E'}" for x in range(-10,141,20)],
    fontsize=6.5)
ax.set_yticks(range(30, 71, 10))
ax.set_yticklabels([f"{y}°N" for y in range(30,71,10)], fontsize=6.5)
ax.tick_params(length=2.5, width=0.35, color="#888888", pad=2)
for sp in ax.spines.values():
    sp.set_linewidth(0.55); sp.set_edgecolor("#555555")

# --- 7. LEYENDA: BARRA LIMPIA DENTRO DEL MAPA ---------------
cmap_leg = LinearSegmentedColormap.from_list("gen",
    [(0.00, (205/255,  18/255,  42/255)),
     (0.35, (228/255,  62/255,  12/255)),
     (0.65, ( 18/255,  55/255, 210/255)),
     (1.00, ( 55/255,  35/255, 180/255))])

# Inset container
ax_in = inset_axes(ax, width="34%", height="22%",
                   loc="lower right", borderpad=0.8)
ax_in.set_facecolor("white")
for sp in ax_in.spines.values():
    sp.set_linewidth(0.4); sp.set_edgecolor("#CCCCCC")
ax_in.tick_params(left=False, bottom=False,
                  labelleft=False, labelbottom=False)
ax_in.set_xlim(0, 1); ax_in.set_ylim(0, 1)

# Barra de color dentro del inset
ax_cb = inset_axes(ax_in, width="76%", height="32%",
                   loc="center", borderpad=0)
cb = ColorbarBase(ax_cb, cmap=cmap_leg,
                  orientation="horizontal",
                  norm=mcolors.Normalize(0, 1))
cb.outline.set_linewidth(0.35)
cb.outline.set_edgecolor("#999999")
ax_cb.tick_params(length=0, labelsize=0)

# Muescas Vindija y Denisova
for pos, lbl, col in [
    (0.35, "Vindija",  (228/255,  62/255,  12/255)),
    (0.65, "Denisova", ( 55/255,  35/255, 180/255)),
]:
    ax_cb.plot([pos], [0.5], "|", transform=ax_cb.transAxes,
               color="white", markersize=6, markeredgewidth=1.0,
               clip_on=False, zorder=5)
    ax_cb.text(pos, -0.6, lbl, transform=ax_cb.transAxes,
               ha="center", va="top", fontsize=5.5,
               color=col, style="italic")

# Etiquetas de polo
ax_cb.text(-0.03, 0.5, "Sima de los Huesos",
           transform=ax_cb.transAxes, ha="right", va="center",
           fontsize=5.8, color=(205/255, 18/255, 42/255),
           style="italic", fontweight="bold")
ax_cb.text(1.03, 0.5, "Ust'-Ishim",
           transform=ax_cb.transAxes, ha="left", va="center",
           fontsize=5.8, color=(18/255, 55/255, 210/255),
           style="italic", fontweight="bold")

# Única línea de leyenda: qué mide el gradiente
ax_cb.text(0.5, 1.9,
           "Nuclear genetic affinity  (IDW · Chr. 1)",
           transform=ax_cb.transAxes, ha="center", va="bottom",
           fontsize=5.5, color="#555555")

# --- 8. BARRA DE ESCALA ------------------------------------
ax.plot([-12,8],[28.7,28.7], color="#222222", lw=1.0,
        zorder=8, solid_capstyle="butt")
ax.plot([-12,-12],[28.3,29.1], color="#222222", lw=0.65, zorder=8)
ax.plot([8,8],    [28.3,29.1], color="#222222", lw=0.65, zorder=8)
ax.text(-2, 27.95, "~2 000 km", ha="center", va="top",
        fontsize=6.0, color="#222222", zorder=8)

# --- 9. EXPORTAR --------------------------------------------
for fmt, d in [("png",600),("pdf",300),("tiff",600)]:
    fname = f"mapa_gradiente_final5.{fmt}"
    kw = dict(dpi=d, bbox_inches="tight", facecolor="white")
    if fmt == "tiff":
        kw["pil_kwargs"] = {"compression": "tiff_lzw"}
    plt.savefig(fname, **kw)
    print(f"  → {fname}")

print("\nDone.")
