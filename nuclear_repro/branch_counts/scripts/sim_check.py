#!/usr/bin/env python3
"""Comprobación del recuento por rama sobre datos simulados con topología conocida.

Simula una estructura neandertal en la que Denisova 17 es hermano de Denisova 5, y
Vindija 33.19 lo es de Chagyrskaya 8 (la topología que implican las cifras de Massilani
et al. 2026). Luego aplica branch_counts con X = D17, A = D5 y B = Vi33.19 y comprueba
que el orden de los tres estadísticos es el esperado por la topología:

  f_AB  alto   (rama neandertal profunda: D17 debe compartir casi todo)
  f_A   medio  (específico de D5: D17 comparte lo del clado D5+D17)
  f_B   bajo   (específico de Vindija: D17 no debería compartir casi nada)

No reproduce los valores de Massilani: la demografía real no se conoce. Comprueba que el
código ordena correctamente los tres conjuntos y que f_A > f_B, que es la firma de la
topología. Un código que confundiera los conjuntos A y B daría el orden invertido.

uso: python3 scripts/sim_check.py [cromosomas]
"""
import sys
from pathlib import Path

import msprime
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from branch_counts import branch_counts  # noqa: E402

GEN = 29.0
CHR_LEN = {"21": 48_129_895, "22": 51_304_566}
MU, REC = 1.25e-8, 1e-8
BLOCK = 5_000_000


def ka(x):
    return x * 1000.0 / GEN


def demography():
    d = msprime.Demography()
    for name, n in (("AFR", 15000), ("DEN", 3000), ("NEA_W", 2500), ("NEA_E", 2500),
                    ("NEA", 3000), ("ARC", 5000), ("HOM", 10000)):
        d.add_population(name=name, initial_size=n)
    # (D5, D17) y (Vi, Chag) como dos subclados neandertales
    d.add_population_split(time=ka(130), derived=["NEA_W", "NEA_E"], ancestral="NEA")
    d.add_population_split(time=ka(470), derived=["DEN", "NEA"], ancestral="ARC")
    d.add_population_split(time=ka(600), derived=["AFR", "ARC"], ancestral="HOM")
    return d


SAMPLES = [
    msprime.SampleSet(1, population="AFR", time=0),          # Mbuti
    msprime.SampleSet(1, population="DEN", time=ka(65)),     # Denisova 3
    msprime.SampleSet(1, population="NEA_E", time=ka(118)),  # Denisova 5
    msprime.SampleSet(1, population="NEA_E", time=ka(110)),  # Denisova 17
    msprime.SampleSet(1, population="NEA_W", time=ka(49)),   # Vindija 33.19
    msprime.SampleSet(1, population="NEA_W", time=ka(77)),   # Chagyrskaya 8
]
NAMES = ["mbuti", "den3", "d5", "d17", "vindija", "chag8"]


def main(chroms):
    d = demography()
    cols, blocks = {n: [] for n in NAMES}, []
    for i, c in enumerate(chroms):
        ts = msprime.sim_ancestry(samples=SAMPLES, demography=d, ploidy=1,
                                  sequence_length=CHR_LEN[c], recombination_rate=REC,
                                  random_seed=910_000 + i)
        ts = msprime.sim_mutations(ts, rate=MU, random_seed=450_000 + i,
                                   model=msprime.BinaryMutationModel())
        g = ts.genotype_matrix()          # 0 = ancestral, 1 = derivado
        pos = ts.sites_position.astype(np.int64)
        for j, n in enumerate(NAMES):
            cols[n].append(g[:, j].astype(np.int8))
        blocks.append(i * 1000 + (pos // BLOCK).astype(np.int64))
        print(f"  cromosoma {c}: {g.shape[0]:,} sitios segregantes")
    arr = {n: np.concatenate(v) for n, v in cols.items()}
    blk = np.concatenate(blocks)
    _, blk = np.unique(blk, return_inverse=True)

    r = branch_counts(arr["d17"], arr["d5"], arr["vindija"], arr["mbuti"], blocks=blk)
    print("\nX = Denisova 17, A = Denisova 5, B = Vindija 33.19")
    print("conjunto                      sitios        f      EE")
    etq = {"AB": "AB  rama profunda", "A": "A   exclusivos de D5", "B": "B   exclusivos de Vindija"}
    for k in ("AB", "A", "B"):
        v = r[k]
        ee = f"{v['ee']:.4f}" if v["ee"] is not None and np.isfinite(v["ee"]) else "n/d"
        print(f"{etq[k]:28s} {v['n']:8,}  {v['f']:.3f}  {ee}")

    ok = []
    ok.append(("f_AB es el mayor de los tres", r["AB"]["f"] > r["A"]["f"] and r["AB"]["f"] > r["B"]["f"]))
    ok.append(("f_A > f_B (firma de la topología)", r["A"]["f"] > r["B"]["f"]))
    ok.append(("los tres conjuntos tienen datos", all(r[k]["n"] > 100 for k in r)))
    # control: con X = Chagyrskaya el orden debe invertirse, porque Chag va con Vindija
    rc = branch_counts(arr["chag8"], arr["d5"], arr["vindija"], arr["mbuti"])
    ok.append(("control X = Chagyrskaya invierte f_A y f_B", rc["B"]["f"] > rc["A"]["f"]))
    print(f"\ncontrol X = Chagyrskaya 8:  f_A = {rc['A']['f']:.3f}   f_B = {rc['B']['f']:.3f}")
    print()
    bad = [n for n, v in ok if not v]
    for n, v in ok:
        print(f"  {'ok  ' if v else 'FALLA'} {n}")
    print("\nRESULTADO: " + ("todas las comprobaciones pasan" if not bad else f"fallan: {bad}"))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:] or ["21", "22"]))
