#!/usr/bin/env python3
"""Paso 6: Figura 1 dibujada directamente desde los árboles de IQ-TREE (sin teclear topologías).

Sustituye a documents/JHE/submission/figure1_mtDNA.R y figures/tanglegramSima.R, que
escribían el árbol a mano como texto Newick (ver CORRECTION.md).

Lee results/trees/iq_{full,coding}_mfp_ufboot.contree, enraíza en Pan_outgro y etiqueta cada
nodo con el soporte UFBoot que figura en el propio fichero. Además anota el soporte de
Sima+Denisova en las 1000 réplicas de cada método, leído de results/CLADE_SUPPORT.tsv.

uso: python3 scripts/06_figure.py     (desde mtDNA_repro/, entorno mtphylo)
salidas: results/figures/figure1_mtDNA_rooted.{pdf,png}
"""
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from Bio import Phylo  # noqa: E402

NAMES = {
    "Sima_430ka": "Sima de los Huesos (Femur XIII, ~430 ka)",
    "Denisova3_": "Denisova 3 (~70 ka)",
    "Neandertal": "Vindija 33.16 (~44 ka)",
    "Feldhofer1": "Feldhofer 1 (~40 ka)",
    "Mezmaiskay": "Mezmaiskaya 1 (~60 ka)",
    "Pan_outgro": "Pan troglodytes (grupo externo)",
}
PAN_DRAWN = 0.004  # longitud dibujada de la rama de Pan (la real se indica en la etiqueta)


def load(path):
    tree = Phylo.read(path, "newick")
    tree.root_with_outgroup("Pan_outgro")
    # La raíz parte la rama Pan-homininos en dos; su longitud total es la distancia real.
    # Se dibuja acortada (PAN_DRAWN, repartida entre ambos lados) y se indica el valor real.
    real = sum(c.branch_length or 0 for c in tree.root.clades)
    for c in tree.root.clades:
        c.branch_length = PAN_DRAWN / 2
    return tree, real


def support_table():
    rows = list(csv.DictReader(open("results/CLADE_SUPPORT.tsv"), delimiter="\t"))
    return {r["analisis"]: r for r in rows}


def main():
    sup = support_table()
    out = Path("results/figures")
    out.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(15, 7.2))
    panels = [("full", "Genoma completo (17.150 col.)"), ("coding", "Sin D-loop (15.409 col.)")]
    for ax, (aln, title) in zip(axes, panels):
        tree, pan_len = load(f"results/trees/iq_{aln}_mfp_ufboot.contree")
        model = open(f"results/trees/iq_{aln}_mfp_ufboot.iqtree").read().split(
            "Best-fit model according to BIC: ")[1].split()[0]

        def label(c):
            if c.name == "Pan_outgro":
                return f"{NAMES[c.name]}  [rama acortada; real {pan_len:.3f} sust./sitio]"
            return NAMES.get(c.name, "")

        def blabel(c):
            if c.is_terminal() or c.confidence is None:
                return None
            return f"{c.confidence:.0f}"

        Phylo.draw(tree, axes=ax, do_show=False, label_func=label, branch_labels=blabel)
        for t in ax.texts:
            if t.get_text().strip().startswith("Sima"):
                t.set_color("#B03030")
                t.set_fontweight("bold")
        ax.set_title(f"{title}\nIQ-TREE 3.0.1, {model}, raíz en Pan · nodos: UFBoot (1000 réplicas)",
                     fontsize=10)
        ax.set_ylabel("")
        ax.set_yticks([])
        ax.set_xlabel("Sustituciones por sitio")
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
        lines = ["Sima+Denisova en las réplicas bootstrap:"]
        for k, r in sup.items():
            if k.startswith(f"{aln}_"):
                method = k.split("_", 1)[1].replace("_", " ").replace("stdboot", "bootstrap est.")
                lines.append(f"  {method:32s} {r['%Sima+Denisova']:>5} %")
        ax.text(0.5, -0.17, "\n".join(lines), transform=ax.transAxes, ha="center", va="top",
                fontsize=7.5, family="monospace",
                bbox=dict(boxstyle="round", fc="#F5F5F5", ec="#CCCCCC"))
    fig.suptitle("Figura 1. Filogenia mitocondrial enraizada: Sima de los Huesos es hermano de Denisova",
                 fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    for ext in ("pdf", "png"):
        fig.savefig(out / f"figure1_mtDNA_rooted.{ext}", dpi=200)
    print("Escrito:", *sorted(str(p) for p in out.glob("figure1_mtDNA_rooted.*")))


if __name__ == "__main__":
    main()
