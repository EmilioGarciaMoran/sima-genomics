#!/usr/bin/env python3
"""Paso 5: resume las pruebas de la ampliación y aplica los criterios del README.

uso: python3 scripts/05_ext_summary.py results > results/EXT_SUMMARY.tsv
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

base = Path(sys.argv[1])
rows = defaultdict(list)
for p in sorted(base.glob("TND*/EXT_REPLICATES.tsv")):
    for r in csv.DictReader(open(p), delimiter="\t"):
        rows[(r["escenario"], r["X"], r["nivel"])].append(r)

ORDER = {"femurXIII": 0, "AT5431": 1, "AT5482": 2, "AT5431+AT5482": 3, "full": 4}
frac = lambda z, f: float(np.mean(f(np.nan_to_num(z, nan=0.0))))  # noqa: E731

print("escenario\tX\tnivel\tn_D_mediana\tD_Z>3\tcub_N+cub_D_mediana\tf_N_mediana\tf_D_mediana\tMeyer_Z>3\t"
      "D_|Z|<3\tMeyer_|Z|<3")
res = {}
for k in sorted(rows, key=lambda k: (k[0], k[1], ORDER[k[2]])):
    rs = rows[k]
    zd = np.array([float(r["Z_D"]) for r in rs]); zm = np.array([float(r["Z_Meyer"]) for r in rs])
    nd = np.median([int(r["n_D"]) for r in rs])
    cov = np.median([int(r["cub_N"]) + int(r["cub_D"]) for r in rs])
    fN = np.nanmedian([float(r["f_N"]) for r in rs]); fD = np.nanmedian([float(r["f_D"]) for r in rs])
    a, b = frac(zd, lambda z: z > 3), frac(zm, lambda z: z > 3)
    c, d = frac(zd, lambda z: np.abs(z) < 3), frac(zm, lambda z: np.abs(z) < 3)
    res[k] = (a, b, c, d)
    print(f"{k[0]}\t{k[1]}\t{k[2]}\t{int(nd)}\t{a:.2f}\t{int(cov)}\t{fN:.3f}\t{fD:.3f}\t{b:.2f}\t{c:.2f}\t{d:.2f}")

print()
print("# Veredicto por escenario y volumen (Sima Z>3 en >=90 % y Mbuti |Z|<3 en >=95 %)")
for esc in sorted({k[0] for k in rows}):
    for lvl in ORDER:
        s, m = res.get((esc, "sima", lvl)), res.get((esc, "mbuti", lvl))
        if not (s and m):
            continue
        vD = "SÍ" if s[0] >= 0.90 and m[2] >= 0.95 else "no"
        vM = "SÍ" if s[1] >= 0.90 and m[3] >= 0.95 else "no"
        print(f"# {esc}  {lvl:14s}  D: {vD} ({s[0]:.0%})   Meyer: {vM} ({s[1]:.0%})")

print()
print("# Calibración: cociente f_N/f_D de Sima con genoma completo (Meyer observó 4,8 en AT-5431 y 5,6 en AT-5482)")
for p in sorted(base.glob("TND*/CALIBRATION.tsv")):
    for r in csv.DictReader(open(p), delimiter="\t"):
        if r["X"] == "sima":
            print(f"# {r['escenario']}: f_N/f_D = {float(r['cociente_fN_fD']):.2f}  "
                  f"(f_ND={float(r['f_ND']):.2f}, f_N={float(r['f_N']):.2f}, f_D={float(r['f_D']):.2f})")
