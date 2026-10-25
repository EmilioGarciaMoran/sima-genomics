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


def ee_binom(v):
    """Error típico binomial de k/n. Fiable cuando los bloques son pocos."""
    n, f = v["n"], v["f"]
    return float(np.sqrt(f * (1 - f) / n)) if n and np.isfinite(f) else np.nan


def show(title, r, extra=""):
    print(f"\n{title}{extra}")
    print(f"  {'conjunto':38s} {'sitios':>9s} {'derivados':>10s} {'f':>8s} "
          f"{'EE binom':>9s} {'EE jackk':>9s}")
    for k in ("AB", "A", "B"):
        v = r[k]
        ee = f"{v['ee']:.4f}" if v["ee"] is not None and np.isfinite(v["ee"]) else "n/d"
        eb = ee_binom(v)
        f = f"{v['f']:.4f}" if np.isfinite(v["f"]) else "n/d"
        print(f"  {ETQ[k]:38s} {v['n']:9,} {v['k']:10,} {f:>8s} "
              f"{eb:9.4f} {ee:>9s}")
    return r


def z_diff(r, use_jackknife):
    """Z de la diferencia f_A - f_B.

    Usa el error típico del jackknife solo si hay bloques suficientes; con pocos
    bloques el jackknife es inestable y el binomial es la opción honesta.
    """
    a, b = r["A"], r["B"]
    if use_jackknife and all(x["ee"] and np.isfinite(x["ee"]) for x in (a, b)):
        sa, sb = a["ee"], b["ee"]
    else:
        sa, sb = ee_binom(a), ee_binom(b)
    if not (np.isfinite(sa) and np.isfinite(sb)) or (sa == 0 and sb == 0):
        return np.nan
    return (a["f"] - b["f"]) / np.sqrt(sa ** 2 + sb ** 2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("npz")
    ap.add_argument("--block", type=int, default=1_000_000,
                    help="tamaño de bloque del jackknife en pb (por defecto 1 Mb)")
    ap.add_argument("--exclude-rm", action="store_true",
                    help="descarta posiciones marcadas RM (repetición enmascarada)")
    ap.add_argument("--map20", choices=["todo", "ausente", "uno"], default="todo",
                    help="todo: sin filtro; ausente: solo donde Map20 NO está anotado; "
                         "uno: solo donde Map20 >= 1")
    a = ap.parse_args()

    z = np.load(a.npz, allow_pickle=False)
    pos, n_all = z["pos"], len(z["pos"])
    keep = np.ones(n_all, bool)
    if "map20" in z.files and a.map20 != "todo":
        m = z["map20"]
        keep &= (m < 0) if a.map20 == "ausente" else (m >= 1)
    if a.exclude_rm:
        if "rm" in z.files:
            keep &= ~z["rm"]
        else:
            print("AVISO: el npz no trae el campo rm; vuelve a generarlo con read_vcfs.py")

    pos = pos[keep]
    col = {n: z[n][keep] for n in ("D5", "Vi3319", "D3", "Chag8", "Mbuti")}
    blocks = (pos // a.block).astype(np.int64)
    _, blocks = np.unique(blocks, return_inverse=True)
    nblk = blocks.max() + 1 if len(blocks) else 0
    use_jk = nblk >= 20

    print(f"Fichero: {a.npz}")
    print(f"  modo: {z['meta_mode'][0]}   semilla: {z['meta_seed'][0]}   "
          f"exclude_rm(lectura): {z['meta_exclude_rm'][0]}   min_dp: {z['meta_min_dp'][0]}")
    print(f"  filtros de este análisis: map20={a.map20}, exclude_rm={a.exclude_rm}")
    print(f"  posiciones: {len(pos):,} de {n_all:,}"
          f"   rango: {pos.min():,} - {pos.max():,}" if len(pos) else "  sin posiciones")
    print(f"  bloques de jackknife: {nblk} (de {a.block/1e6:.1f} Mb)")
    if not use_jk:
        print("  AVISO: menos de 20 bloques. El jackknife no es fiable aquí;")
        print("         el veredicto usa el error típico binomial.")

    r_chag = show("PRUEBA PRINCIPAL  X = Chagyrskaya 8, A = Denisova 5, B = Vindija 33.19",
                  branch_counts(col["Chag8"], col["D5"], col["Vi3319"], col["Mbuti"],
                                blocks=blocks))
    r_d3 = show("CONTROL           X = Denisova 3 (denisovano)",
                branch_counts(col["D3"], col["D5"], col["Vi3319"], col["Mbuti"],
                              blocks=blocks))

    print("\n" + "=" * 64)
    print("VEREDICTO")
    print("=" * 64)
    zc = z_diff(r_chag, use_jk)
    ok1 = np.isfinite(r_chag["B"]["f"]) and r_chag["B"]["f"] > r_chag["A"]["f"]
    print(f"  1. Topología: f_B > f_A con Chagyrskaya 8")
    print(f"     f_A = {r_chag['A']['f']:.4f}   f_B = {r_chag['B']['f']:.4f}   "
          f"Z(f_A - f_B) = {zc:.1f}")
    print(f"     -> {'PASA' if ok1 else 'FALLA'}"
          f"{'' if abs(zc) >= 3 or not ok1 else '  (pero la diferencia no es significativa)'}")

    ok2 = r_chag["AB"]["f"] > max(r_chag["A"]["f"], r_chag["B"]["f"])
    print(f"  2. f_AB es el mayor de los tres: {'PASA' if ok2 else 'FALLA'}")

    zd = z_diff(r_d3, use_jk)
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
