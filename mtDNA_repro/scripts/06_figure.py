#!/usr/bin/env python3
"""Paso 6: Figura 1 dibujada directamente desde los árboles de IQ-TREE (sin teclear topologías).

Sustituye a documents/JHE/submission/figure1_mtDNA.R y figures/tanglegramSima.R, que
escribían el árbol a mano como texto Newick (ver CORRECTION.md).

Lee results/trees/iq_{full,coding}_mfp_ufboot.contree, enraíza en Pan_outgro y etiqueta cada
nodo con el soporte UFBoot que figura en el propio fichero. Además anota el soporte de
Sima+Denisova en las 1000 réplicas de cada método, leído de results/CLADE_SUPPORT.tsv.

uso: python3 scripts/06_figure.py [es|en]     (desde mtDNA_repro/, entorno mtphylo; por defecto es)
salidas: results/figures/figure1_mtDNA_rooted.{pdf,png} (es) o figure1_mtDNA_rooted_en.{pdf,png} (en)
"""
import csv
import sys
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
}
TXT = {
    "es": dict(pan="Pan troglodytes (grupo externo)", short="rama acortada; real {:.3f} sust./sitio",
               full="Genoma completo (17.150 col.)", coding="Sin D-loop (15.409 col.)",
               sub="IQ-TREE 3.0.1, {}, raíz en Pan · nodos: UFBoot (1000 réplicas)",
               x="Sustituciones por sitio", box="Sima+Denisova en las réplicas bootstrap:",
               std="bootstrap est.",
               title="Figura 1. Filogenia mitocondrial enraizada: Sima de los Huesos es hermano de Denisova",
               suffix=""),
    "en": dict(pan="Pan troglodytes (outgroup)", short="branch shortened; real {:.3f} subst./site",
               full="Complete genome (17,150 col.)", coding="D-loop removed (15,409 col.)",
               sub="IQ-TREE 3.0.1, {}, rooted on Pan · nodes: UFBoot (1,000 replicates)",
               x="Substitutions per site", box="Sima+Denisova in bootstrap replicates:",
               std="standard bootstrap",
               title="Figure 1. Rooted mitochondrial phylogeny: Sima de los Huesos is sister to Denisova",
               suffix="_en"),
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


def main(lang="es"):
    L = TXT[lang]
    sup = support_table()
    out = Path("results/figures")
    out.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(15, 7.2))
    panels = [("full", L["full"]), ("coding", L["coding"])]
    for ax, (aln, title) in zip(axes, panels):
        tree, pan_len = load(f"results/trees/iq_{aln}_mfp_ufboot.contree")
        model = open(f"results/trees/iq_{aln}_mfp_ufboot.iqtree").read().split(
            "Best-fit model according to BIC: ")[1].split()[0]

        def label(c):
            if c.name == "Pan_outgro":
                return f"{L['pan']}  [{L['short'].format(pan_len)}]"
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
        ax.set_title(f"{title}\n{L['sub'].format(model)}", fontsize=10)
        ax.set_ylabel("")
        ax.set_yticks([])
        ax.set_xlabel(L["x"])
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
        lines = [L["box"]]
        for k, r in sup.items():
            if k.startswith(f"{aln}_"):
                method = k.split("_", 1)[1].replace("_", " ").replace("stdboot", L["std"])
                lines.append(f"  {method:32s} {r['%Sima+Denisova']:>5} %")
        ax.text(0.5, -0.17, "\n".join(lines), transform=ax.transAxes, ha="center", va="top",
                fontsize=7.5, family="monospace",
                bbox=dict(boxstyle="round", fc="#F5F5F5", ec="#CCCCCC"))
    fig.suptitle(L["title"], fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    for ext in ("pdf", "png"):
        fig.savefig(out / f"figure1_mtDNA_rooted{L['suffix']}.{ext}", dpi=200)
    print("Escrito:", *sorted(str(p) for p in out.glob(f"figure1_mtDNA_rooted{L['suffix']}.*")))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "es")
