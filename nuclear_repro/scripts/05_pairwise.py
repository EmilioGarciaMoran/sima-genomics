#!/usr/bin/env python3
"""Paso 5: compone la matriz sitios × especímenes a partir de los vectores y calcula
distancias por pares y estadísticos de alelos derivados.

Entradas (versionadas, generadas en la nube por cloud/run_vectors.sh):
    pairwise/panel.tsv.gz, pairwise/vectors/<espécimen>.vec.gz, pairwise/specimens.tsv

Salidas en results/pairwise/:
  specimens_summary.tsv  sitios con alelo por espécimen
  distances.tsv          distancia p por pares (diferencias / sitios con alelo en ambos), EE por jackknife
  distance_matrix.tsv    la misma, en matriz cuadrada
  derived_sharing.tsv    para cada target T y comparador X del grupo interno (los del grupo del
                         externo, p. ej. otros simios, se excluyen): de los sitios donde X lleva el alelo
                         derivado (≠ grupo externo), fracción en que T también lo lleva
  d_stats.tsv            D(T; X, Y) = (nX - nY) / (nX + nY), con nX = sitios donde T comparte el
                         derivado con X y no con Y, nY al revés. Los errores de T se reparten por
                         igual entre nX y nY, así que D es robusto a la baja calidad de T; |Z| > 3
                         indica que T está más cerca de X (D > 0) o de Y (D < 0)

Cuidado al leer distances.tsv: el panel está seleccionado (solo sitios variables en los
comparadores), así que las distancias son relativas al panel, no divergencias por pb; y
la distancia de un espécimen de baja cobertura incluye su tasa de error. La fila de cada
target con el grupo externo mide ese exceso (comparar con la de un comparador al mismo simio).

    python3 scripts/05_pairwise.py [pairwise] [results]
"""
from __future__ import annotations

import csv
import itertools
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import die, load_config  # noqa: E402
from pairwise_io import Panel, read_vector  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent


def read_specimens(path: Path):
    lines = [l for l in path.read_text().splitlines() if l.strip() and not l.startswith("##")]
    return list(csv.DictReader(lines, delimiter="\t"))


def block_ids(panel: Panel, size: int) -> np.ndarray:
    ids, off = [], 0
    for pos in panel.chroms.values():
        b = (pos // size).astype(np.int64)
        ids.append(b + off)
        off += int(b.max()) + 1 if len(b) else 0
    return np.unique(np.concatenate(ids), return_inverse=True)[1] if ids else np.zeros(0, np.int64)


def per_block(mask: np.ndarray, blocks: np.ndarray, nb: int) -> np.ndarray:
    return np.bincount(blocks[mask], minlength=nb).astype(float)


def jackknife(stat, *by_block):
    """Jackknife delete-one-block (bloques sin datos excluidos). Devuelve (estimación, EE)."""
    tot = [b.sum() for b in by_block]
    est = stat(*tot)
    used = np.nonzero(sum(by_block) > 0)[0]
    g = len(used)
    if g < 2 or not np.isfinite(est):
        return est, float("nan")
    reps = np.array([stat(*(t - b[i] for t, b in zip(tot, by_block))) for i in used])
    reps = reps[np.isfinite(reps)]
    return est, float(np.sqrt((g - 1) / g * np.sum((reps - reps.mean()) ** 2)))


def ratio(a, b):
    return a / b if b else float("nan")


def dstat(a, b):
    return (a - b) / (a + b) if (a + b) else float("nan")


def main(pw: Path, out: Path) -> None:
    cfg = load_config(ROOT / "config.env")
    size = int(cfg["BLOCK_SIZE"])
    panel = Panel.read(pw / "panel.tsv.gz")
    specs = read_specimens(pw / "specimens.tsv")
    names, groups, roles, cols = [], {}, {}, []
    for s in specs:
        vec = pw / "vectors" / f"{s['name']}.vec.gz"
        if not vec.exists():
            print(f"  aviso: falta el vector de {s['name']}; se omite")
            continue
        try:
            _, v = read_vector(vec, panel)
        except ValueError as e:
            die(str(e))
        names.append(s["name"])
        groups[s["name"]], roles[s["name"]] = s["group"], s["role"]
        cols.append(v)
    if len(names) < 2:
        die("hacen falta al menos dos vectores")
    M = np.stack(cols, axis=1)           # sitios × especímenes; 0 = N, 1 = REF, 2 = ALT
    called = M > 0
    blocks = block_ids(panel, size)
    nb = int(blocks.max()) + 1
    idx = {n: i for i, n in enumerate(names)}
    out.mkdir(parents=True, exist_ok=True)
    print(f"Matriz: {panel.n:,} sitios × {len(names)} especímenes, {nb} bloques de {size:,} pb")

    with open(out / "specimens_summary.tsv", "w") as fh:
        fh.write("specimen\tgroup\trole\tsites_called\tfrac_called\n")
        for n in names:
            c = int(called[:, idx[n]].sum())
            fh.write(f"{n}\t{groups[n]}\t{roles[n]}\t{c}\t{c / panel.n:.4f}\n")

    D = np.full((len(names), len(names)), np.nan)
    with open(out / "distances.tsv", "w") as fh:
        fh.write("specimen_1\tspecimen_2\tsites_compared\tdifferences\tdistance\tse_jackknife\n")
        for a, b in itertools.combinations(names, 2):
            i, j = idx[a], idx[b]
            both = called[:, i] & called[:, j]
            diff = both & (M[:, i] != M[:, j])
            d, se = jackknife(ratio, per_block(diff, blocks, nb), per_block(both, blocks, nb))
            D[i, j] = D[j, i] = d
            fh.write(f"{a}\t{b}\t{int(both.sum())}\t{int(diff.sum())}\t{d:.6f}\t{se:.6f}\n")
    np.fill_diagonal(D, 0.0)
    with open(out / "distance_matrix.tsv", "w") as fh:
        fh.write("\t" + "\t".join(names) + "\n")
        for n in names:
            fh.write(n + "\t" + "\t".join(f"{x:.6f}" for x in D[idx[n]]) + "\n")

    outgroups = [n for n in names if roles[n] == "outgroup"]
    targets = [n for n in names if roles[n] == "target"]
    comps = [n for n in names if roles[n] == "comparator"]
    if not outgroups:
        print("  sin grupo externo: se omiten alelos derivados y D")
        return
    o = idx[outgroups[0]]
    O = M[:, o]
    # Polarizar con el grupo externo solo tiene sentido para comparadores del grupo interno:
    # se excluyen los del mismo grupo que el externo (otros simios)
    ingroup = [c for c in comps if groups[c] != groups[outgroups[0]]]
    with open(out / "derived_sharing.tsv", "w") as fh:
        fh.write("target\tcomparator\toutgroup\tsites_derived_in_comparator\tshared\tfraction\tse_jackknife\n")
        for t in targets:
            for x in ingroup:
                T, X = M[:, idx[t]], M[:, idx[x]]
                base = called[:, o] & called[:, idx[x]] & called[:, idx[t]] & (X != O)
                shared = base & (T == X)
                f, se = jackknife(ratio, per_block(shared, blocks, nb), per_block(base, blocks, nb))
                fh.write(f"{t}\t{x}\t{outgroups[0]}\t{int(base.sum())}\t{int(shared.sum())}\t{f:.6f}\t{se:.6f}\n")
    with open(out / "d_stats.tsv", "w") as fh:
        fh.write("target\tX\tY\toutgroup\tnX\tnY\tD\tse_jackknife\tZ\n")
        for t in targets:
            for x, y in itertools.combinations(ingroup, 2):
                T, X, Y = M[:, idx[t]], M[:, idx[x]], M[:, idx[y]]
                ok = called[:, o] & called[:, idx[x]] & called[:, idx[y]] & called[:, idx[t]]
                nx = ok & (T == X) & (X != O) & (Y == O)
                ny = ok & (T == Y) & (Y != O) & (X == O)
                d, se = jackknife(dstat, per_block(nx, blocks, nb), per_block(ny, blocks, nb))
                z = d / se if se and np.isfinite(se) and se > 0 else float("nan")
                fh.write(f"{t}\t{x}\t{y}\t{outgroups[0]}\t{int(nx.sum())}\t{int(ny.sum())}\t"
                         f"{d:.4f}\t{se:.4f}\t{z:.2f}\n")
    for t in targets:
        close = sorted((D[idx[t], idx[c]], c) for c in comps + outgroups)
        print(f"  {t}: más cercano a " + ", ".join(f"{c} ({d:.4f})" for d, c in close[:3]))


if __name__ == "__main__":
    a = sys.argv[1:] + [None] * 2
    main(Path(a[0] or ROOT / "pairwise"), Path(a[1] or ROOT / "results" / "pairwise"))
