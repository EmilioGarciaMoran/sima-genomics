#!/usr/bin/env python3
"""Aplica el recuento por rama a los alelos leídos y ejecuta las pruebas de validación.

Prueba principal (topología, preregistrada en README.md):
    X = Chagyrskaya 8, A = Denisova 5, B = Vindija 33.19
    La topología establecida exige f_B > f_A: Chag8 está más emparentada con Vindija
    que con el neandertal del Altái.

Control:
    X = Denisova 3 (denisovano, grupo externo a los dos neandertales)
    No debe inclinarse hacia ninguno: f_A y f_B parecidos, y ambos mucho menores que f_AB.

uso:
    python3 scripts/run_validation.py results/chr21_alleles.npz [--block 1000000]
"""
import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from branch_counts import branch_counts  # noqa: E402

ETQ = {"AB": "AB  rama profunda (D5 y Vi derivados)",
       "A": "A   exclusivos de D5",
       "B": "B   exclusivos de Vindija"}


def show(title, r, extra=""):
    print(f"\n{title}{extra}")
    print(f"  {'conjunto':38s} {'sitios':>9s} {'derivados':>10s} {'f':>8s} {'EE':>8s}")
    for k in ("AB", "A", "B"):
        v = r[k]
        ee = f"{v['ee']:.4f}" if v["ee"] is not None and np.isfinite(v["ee"]) else "n/d"
        f = f"{v['f']:.4f}" if np.isfinite(v["f"]) else "n/d"
        print(f"  {ETQ[k]:38s} {v['n']:9,} {v['k']:10,} {f:>8s} {ee:>8s}")
    return r


def z_diff(r):
    """Z de la diferencia f_A - f_B, con los EE del jackknife."""
    a, b = r["A"], r["B"]
    if not (a["ee"] and b["ee"] and np.isfinite(a["ee"]) and np.isfinite(b["ee"])):
        return np.nan
    return (a["f"] - b["f"]) / np.sqrt(a["ee"] ** 2 + b["ee"] ** 2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("npz")
    ap.add_argument("--block", type=int, default=1_000_000,
                    help="tamaño de bloque del jackknife en pb (por defecto 1 Mb)")
    a = ap.parse_args()

    z = np.load(a.npz, allow_pickle=False)
    pos = z["pos"]
    blocks = (pos // a.block).astype(np.int64)
    _, blocks = np.unique(blocks, return_inverse=True)
    nblk = blocks.max() + 1

    print(f"Fichero: {a.npz}")
    print(f"  modo: {z['meta_mode'][0]}   semilla: {z['meta_seed'][0]}   "
          f"exclude_rm: {z['meta_exclude_rm'][0]}   min_dp: {z['meta_min_dp'][0]}")
    print(f"  posiciones: {len(pos):,}   rango: {pos.min():,} - {pos.max():,}")
    print(f"  bloques de jackknife: {nblk} (de {a.block/1e6:.1f} Mb)")
    if nblk < 20:
        print("  AVISO: menos de 20 bloques; el error típico del jackknife será poco fiable.")

    r_chag = show("PRUEBA PRINCIPAL  X = Chagyrskaya 8, A = Denisova 5, B = Vindija 33.19",
                  branch_counts(z["Chag8"], z["D5"], z["Vi3319"], z["Mbuti"], blocks=blocks))
    r_d3 = show("CONTROL           X = Denisova 3 (denisovano)",
                branch_counts(z["D3"], z["D5"], z["Vi3319"], z["Mbuti"], blocks=blocks))

    print("\n" + "=" * 64)
    print("VEREDICTO")
    print("=" * 64)
    zc = z_diff(r_chag)
    ok1 = np.isfinite(r_chag["B"]["f"]) and r_chag["B"]["f"] > r_chag["A"]["f"]
    print(f"  1. Topología: f_B > f_A con Chagyrskaya 8")
    print(f"     f_A = {r_chag['A']['f']:.4f}   f_B = {r_chag['B']['f']:.4f}   "
          f"Z(f_A - f_B) = {zc:.1f}")
    print(f"     -> {'PASA' if ok1 else 'FALLA'}"
          f"{'' if abs(zc) >= 3 or not ok1 else '  (pero la diferencia no es significativa)'}")

    ok2 = r_chag["AB"]["f"] > max(r_chag["A"]["f"], r_chag["B"]["f"])
    print(f"  2. f_AB es el mayor de los tres: {'PASA' if ok2 else 'FALLA'}")

    zd = z_diff(r_d3)
    ok3 = np.isfinite(zd) and abs(zd) < 3
    print(f"  3. Control denisovano sin inclinación: Z(f_A - f_B) = {zd:.1f}")
    print(f"     -> {'PASA' if ok3 else 'FALLA (D3 se inclina hacia un neandertal)'}")

    nmin = min(r_chag[k]["n"] for k in ("AB", "A", "B"))
    print(f"\n  Conjunto más pequeño: {nmin:,} sitios")
    for n, f in ((r_chag["AB"]["n"], 0.91), (r_chag["A"]["n"], 0.45), (r_chag["B"]["n"], 0.20)):
        if n:
            print(f"    con n={n:,}, IC95 sobre una fracción de {f:.0%} es "
                  f"+/- {1.96*np.sqrt(f*(1-f)/n)*100:.2f} puntos")
    print("    (la tolerancia preregistrada es +/- 3 puntos)")
    print()
    print("TODAS LAS PRUEBAS PASAN" if (ok1 and ok2 and ok3) else "HAY PRUEBAS QUE FALLAN")


if __name__ == "__main__":
    main()
