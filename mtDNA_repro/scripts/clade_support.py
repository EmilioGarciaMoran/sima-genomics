#!/usr/bin/env python3
"""Cuenta, en un fichero de árboles de réplicas bootstrap (Newick, uno por línea),
la frecuencia de las tres agrupaciones posibles de Sima respecto a Denisova y los Neandertales,
con Pan como grupo externo. La frecuencia de una bipartición no depende de dónde se dibuje la raíz.

uso: clade_support.py <arboles_bootstrap> <etiqueta>   -> imprime una línea TSV
"""
import re
import sys

NEAND = {"Neandertal", "Feldhofer1", "Mezmaiskay"}
ALL = {"Sima_430ka", "Denisova3_", "Pan_outgro"} | NEAND
CLADES = {
    "Sima+Denisova": {"Sima_430ka", "Denisova3_"},
    "Sima+Neandertales": {"Sima_430ka"} | NEAND,
    "Denisova+Neandertales": {"Denisova3_"} | NEAND,
}


def bipartitions(newick):
    """Conjuntos de hojas bajo cada nodo interno (sin longitudes ni soportes)."""
    s = re.sub(r":[0-9.eE+-]+", "", newick.strip().rstrip(";"))
    s = re.sub(r"\)[0-9.]+", ")", s)
    stack, splits, token = [], [], ""
    for ch in s:
        if ch == "(":
            stack.append(set())
        elif ch in ",)":
            if token:
                stack[-1].add(token.strip())
                token = ""
            if ch == ")":
                node = stack.pop()
                splits.append(frozenset(node))
                if stack:
                    stack[-1] |= node
        else:
            token += ch
    out = set()
    for sp in splits:
        if 1 < len(sp) < len(ALL) - 1:
            # normalizar: lado que no contiene a Pan (equivale a enraizar en Pan)
            out.add(sp if "Pan_outgro" not in sp else frozenset(ALL - sp))
    return out


def main(path, label):
    trees = [l for l in open(path) if l.strip()]
    counts = {k: 0 for k in CLADES}
    for t in trees:
        bp = bipartitions(t)
        for k, v in CLADES.items():
            counts[k] += frozenset(v) in bp
    n = len(trees)
    print("\t".join([label, str(n)] + [f"{100 * counts[k] / n:.1f}" for k in CLADES]))


if __name__ == "__main__":
    if sys.argv[1] == "--header":
        print("\t".join(["analisis", "n_replicas"] + [f"%{k}" for k in CLADES]))
    else:
        main(sys.argv[1], sys.argv[2])
