#!/usr/bin/env python3
"""Une los npz por cromosoma del recorrido genómico en uno solo y resume los conjuntos.

    python3 scripts/merge_genome.py results/genome --out results/genoma_inf.npz
"""
import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from branch_counts import branch_counts  # noqa: E402

COLS = ["D5", "Vi3319", "D3", "Chag8", "Mbuti"]

ap = argparse.ArgumentParser()
ap.add_argument("dir")
ap.add_argument("--out", default="results/genoma_inf.npz")
ap.add_argument("--block", type=int, default=5_000_000)
a = ap.parse_args()

files = sorted(Path(a.dir).glob("chr*_inf.npz"),
               key=lambda p: int(p.name.split("_")[0][3:]))
if not files:
    sys.exit(f"ERROR: no hay chr*_inf.npz en {a.dir}")

acc = {c: [] for c in COLS}
pos, chrom, canc, der = [], [], [], []
print(f"{'cr':>3s} {'sitios':>12s}")
for f in files:
    z = np.load(f, allow_pickle=False)
    c = int(str(z["chrom"][0]))
    n = len(z["pos"])
    print(f"{c:3d} {n:12,}")
    pos.append(z["pos"]); canc.append(z["canc"]); der.append(z["der"])
    chrom.append(np.full(n, c, dtype=np.int8))
    for k in COLS:
        acc[k].append(z[k])

out = {"pos": np.concatenate(pos), "chrom": np.concatenate(chrom),
       "canc": np.concatenate(canc), "der": np.concatenate(der)}
for k in COLS:
    out[k] = np.concatenate(acc[k])
Path(a.out).parent.mkdir(parents=True, exist_ok=True)
np.savez_compressed(a.out, **out)
N = len(out["pos"])
print(f"\nTotal: {N:,} sitios informativos en {len(files)} cromosomas -> {a.out}")

blocks = out["chrom"].astype(np.int64) * 1000 + (out["pos"] // a.block)
_, blocks = np.unique(blocks, return_inverse=True)
print(f"Bloques de jackknife: {blocks.max()+1} (de {a.block/1e6:.0f} Mb)")

r = branch_counts(out["Chag8"], out["D5"], out["Vi3319"], out["Mbuti"], blocks=blocks)
print("\nConjuntos de referencia (genoma completo):")
for k, etq in (("AB", "rama profunda D5+Vi"), ("A", "exclusivos de D5"),
               ("B", "exclusivos de Vindija")):
    v = r[k]
    print(f"  {k:3s} {etq:24s} {v['n']:10,} sitios")
print("\nEste es el universo del que saldrán los sitios de Sima: de cada conjunto,")
print("Sima cubrirá la fracción correspondiente a sus bases alineadas.")
