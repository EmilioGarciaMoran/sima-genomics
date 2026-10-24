#!/usr/bin/env python3
"""Paso 3: resume las réplicas y aplica el criterio de decisión fijado en el README.

uso: python3 scripts/03_summary.py results/*/D_REPLICATES.tsv > results/SUMMARY.tsv
Criterio (escenario principal T_ND = 470 ka): apropiada si Sima da Z > 3 en >= 90 % de réplicas
y Mbuti da |Z| < 3 en >= 95 %.
"""
import csv
import sys
from collections import defaultdict

import numpy as np

rows = defaultdict(list)
for path in sys.argv[1:]:
    for r in csv.DictReader(open(path), delimiter="\t"):
        rows[(r["escenario"], r["X"], r["nivel"])].append(r)

ORDER = {"femurXIII": 0, "AT5431": 1, "AT5482": 2, "full": 3}
print("escenario\tX\tnivel\tn_inf_mediana\tD_mediana\tD_p2.5\tD_p97.5\tEE_mediana\t"
      "frac_Z>3\tfrac_|Z|<3\tveredicto")
summary = {}
for (esc, x, lvl) in sorted(rows, key=lambda k: (k[0], k[1], ORDER[k[2]])):
    rs = rows[(esc, x, lvl)]
    D = np.array([float(r["D"]) for r in rs])
    Z = np.array([float(r["Z"]) for r in rs])
    EE = np.array([float(r["EE"]) for r in rs])
    n = np.array([int(r["n_informativos"]) for r in rs])
    zpos = np.mean(np.nan_to_num(Z, nan=0) > 3)
    zns = np.mean(np.abs(np.nan_to_num(Z, nan=0)) < 3)
    summary[(esc, x, lvl)] = (zpos, zns)
    v = ""
    if x == "sima":
        v = "potencia suficiente" if zpos >= 0.90 else "potencia INSUFICIENTE"
    elif x == "mbuti":
        v = "sin falsos positivos" if zns >= 0.95 else "FALSOS POSITIVOS"
    elif x == "chag":
        v = "control positivo OK" if zpos >= 0.90 else "control positivo FALLA"
    print(f"{esc}\t{x}\t{lvl}\t{int(np.median(n))}\t{np.nanmedian(D):.4f}\t"
          f"{np.nanpercentile(D, 2.5):.4f}\t{np.nanpercentile(D, 97.5):.4f}\t{np.nanmedian(EE):.4f}\t"
          f"{zpos:.2f}\t{zns:.2f}\t{v}")

print()
print("# Decisión (escenario principal T_ND=470ka)")
for lvl in ["femurXIII", "AT5431", "AT5482", "full"]:
    s = summary.get(("TND470", "sima", lvl)); m = summary.get(("TND470", "mbuti", lvl))
    if s and m:
        ok = s[0] >= 0.90 and m[1] >= 0.95
        print(f"# {lvl}: {'D-statistics APROPIADAS' if ok else 'D-statistics NO apropiadas'} "
              f"(Sima Z>3 en {s[0]:.0%}; Mbuti |Z|<3 en {m[1]:.0%})")
