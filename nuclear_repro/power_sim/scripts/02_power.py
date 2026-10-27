#!/usr/bin/env python3
"""Paso 2: D-statistics con submuestreo "tipo Sima" sobre los genomas simulados.

uso: python3 scripts/02_power.py <dir_simulacion> <etiqueta_escenario> <dir_resultados> [replicas]

D(H1 = Denisova 3, H2 = Vindija; X, O = ancestral) = (ABBA - BABA) / (ABBA + BABA)
D > 0  <=>  X comparte más alelos derivados con Vindija (afinidad neandertal).
Todas las muestras pseudohaploides (un alelo al azar por sitio). Solo X se submuestrea.
Error típico: jackknife de bloques ponderado (Busing et al. 1999), bloques de 5 Mb.
"""
import json
import sys
from pathlib import Path

import numpy as np

HAP = {"mbuti": (0, 1), "den3": (2, 3), "vindija": (4, 5), "chag": (6, 7), "sima": (8, 9)}
X_SET = ["sima", "chag", "mbuti"]
LEVELS = {"femurXIII": 225_329, "AT5431": 1_186_442, "AT5482": 2_015_167, "full": None}
BLOCK = 5_000_000


def load(simdir):
    pos, chrom, hap = [], [], []
    for i in range(1, 23):
        z = np.load(Path(simdir) / f"chr{i}.npz")
        pos.append(z["pos"]); hap.append(z["hap"]); chrom.append(np.full(len(z["pos"]), i, np.int16))
    meta = json.load(open(Path(simdir) / "meta.json"))
    return np.concatenate(pos), np.concatenate(chrom), np.concatenate(hap), meta


def pseudohaploid(hap, who, rng):
    a, b = HAP[who]
    pick = rng.random(hap.shape[0]) < 0.5
    return np.where(pick, hap[:, a], hap[:, b])


def jackknife(abba_b, baba_b):
    """D y error típico por jackknife de bloques ponderado (Busing et al. 1999)."""
    m = abba_b + baba_b
    keep = m > 0
    A, B, m = abba_b[keep].astype(float), baba_b[keep].astype(float), m[keep].astype(float)
    n, g = m.sum(), len(m)
    if n == 0:
        return np.nan, np.nan, 0
    d = (A.sum() - B.sum()) / n
    if g < 2:
        return d, np.nan, int(n)
    d_minus = ((A.sum() - A) - (B.sum() - B)) / (n - m)
    h = n / m
    d_j = g * d - np.sum((1 - m / n) * d_minus)
    tau = h * d - (h - 1) * d_minus
    var = np.mean((tau - d_j) ** 2 / (h - 1))
    return d, float(np.sqrt(var)), int(n)


def main(simdir, label, outdir, nrep=100):
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    pos, chrom, hap, meta = load(simdir)
    G = meta["genome_length"]
    seed = int(meta["seed"]) * 1000 + 7
    rng = np.random.default_rng(seed)
    h1 = pseudohaploid(hap, "den3", rng)
    h2 = pseudohaploid(hap, "vindija", rng)
    S = h1 != h2                                   # solo estos sitios pueden ser ABBA o BABA
    blk = (chrom.astype(np.int64) * 10_000 + (pos // BLOCK).astype(np.int64))[S]
    _, blk = np.unique(blk, return_inverse=True)
    nb = blk.max() + 1
    h1S, hapS = h1[S], hap[S]
    abba_pat = (h1S == 0)                           # H1 ancestral, H2 derivado
    rows = []
    for x in X_SET:
        a, b = HAP[x]
        for lvl, bases in LEVELS.items():
            p = 1.0 if bases is None else bases / G
            for r in range(nrep):
                xa = np.where(rng.random(len(h1S)) < 0.5, hapS[:, a], hapS[:, b])
                kept = rng.random(len(h1S)) < p if p < 1 else np.ones(len(h1S), bool)
                der = kept & (xa == 1)
                abba_b = np.bincount(blk[der & abba_pat], minlength=nb)
                baba_b = np.bincount(blk[der & ~abba_pat], minlength=nb)
                d, se, n = jackknife(abba_b, baba_b)
                z = d / se if se and se > 0 else np.nan
                rows.append((label, x, lvl, r, int(abba_b.sum()), int(baba_b.sum()), n, d, se, z))
    with open(out / "D_REPLICATES.tsv", "w") as f:
        f.write("escenario\tX\tnivel\treplica\tABBA\tBABA\tn_informativos\tD\tEE\tZ\n")
        for t in rows:
            f.write("\t".join(f"{v:.5g}" if isinstance(v, float) else str(v) for v in t) + "\n")
    print(f"{label}: {len(rows)} réplicas; sitios con H1≠H2: {int(S.sum()):,}; bloques: {nb}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]) if len(sys.argv) > 4 else 100)
