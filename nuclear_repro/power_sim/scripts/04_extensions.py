#!/usr/bin/env python3
"""Paso 4: pruebas de la ampliación (README): A) volumen combinado AT5431+AT5482 con la D;
B) estadístico de Meyer 2016 (alelos derivados por rama); C) calibración del escenario.

uso: python3 scripts/04_extensions.py <dir_simulacion> <etiqueta> <dir_resultados> [replicas]
salidas: <dir_resultados>/EXT_REPLICATES.tsv y <dir_resultados>/CALIBRATION.tsv
"""
import importlib.util
import sys
from pathlib import Path

import numpy as np

_spec = importlib.util.spec_from_file_location("p2", Path(__file__).with_name("02_power.py"))
p2 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(p2)

LEVELS = {"femurXIII": 225_329, "AT5431": 1_186_442, "AT5482": 2_015_167,
          "AT5431+AT5482": 3_201_609, "full": None}
X_SET = ["sima", "chag", "mbuti"]


def two_prop_z(k1, n1, k2, n2):
    if n1 == 0 or n2 == 0:
        return np.nan
    p1, p2_, p = k1 / n1, k2 / n2, (k1 + k2) / (n1 + n2)
    se = np.sqrt(p * (1 - p) * (1 / n1 + 1 / n2))
    return (p1 - p2_) / se if se > 0 else np.nan


def main(simdir, label, outdir, nrep=100):
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    pos, chrom, hap, meta = p2.load(simdir)
    G = meta["genome_length"]
    rng = np.random.default_rng(int(meta["seed"]) * 1000 + 11)
    den = p2.pseudohaploid(hap, "den3", rng)
    vin = p2.pseudohaploid(hap, "vindija", rng)
    mbu = p2.pseudohaploid(hap, "mbuti", rng)
    # Ramas de Meyer (estado ancestral = 0)
    br = {"ND": (den == 1) & (vin == 1) & (mbu == 0),
          "N": (vin == 1) & (den == 0) & (mbu == 0),
          "D": (den == 1) & (vin == 0) & (mbu == 0)}
    # D-statistic (sitios con H1 != H2), igual que en el paso 2
    S = den != vin
    blk = (chrom.astype(np.int64) * 10_000 + (pos // p2.BLOCK).astype(np.int64))
    _, blk_all = np.unique(blk, return_inverse=True)
    nb = blk_all.max() + 1
    rows, calib = [], []
    for x in X_SET:
        a, b = p2.HAP[x]
        for lvl, bases in LEVELS.items():
            p = 1.0 if bases is None else bases / G
            for r in range(nrep):
                xa = np.where(rng.random(len(pos)) < 0.5, hap[:, a], hap[:, b])
                kept = rng.random(len(pos)) < p if p < 1 else np.ones(len(pos), bool)
                der = kept & (xa == 1)
                # A) D
                abba_b = np.bincount(blk_all[der & S & (den == 0)], minlength=nb)
                baba_b = np.bincount(blk_all[der & S & (den == 1)], minlength=nb)
                d, se, n = p2.jackknife(abba_b, baba_b)
                zd = d / se if se and se > 0 else np.nan
                # B) Meyer: fracción derivada por rama entre los sitios cubiertos
                cov = {k: int((kept & m).sum()) for k, m in br.items()}
                dk = {k: int((der & m).sum()) for k, m in br.items()}
                f = {k: (dk[k] / cov[k] if cov[k] else np.nan) for k in br}
                zm = two_prop_z(dk["N"], cov["N"], dk["D"], cov["D"])
                rows.append((label, x, lvl, r, n, d, zd, cov["ND"], cov["N"], cov["D"],
                             f["ND"], f["N"], f["D"], zm))
                if lvl == "full" and r == 0:
                    calib.append((label, x, f["ND"], f["N"], f["D"], f["N"] / f["D"] if f["D"] else np.nan))
    with open(out / "EXT_REPLICATES.tsv", "w") as fh:
        fh.write("escenario\tX\tnivel\treplica\tn_D\tD\tZ_D\tcub_ND\tcub_N\tcub_D\t"
                 "f_ND\tf_N\tf_D\tZ_Meyer\n")
        for t in rows:
            fh.write("\t".join(f"{v:.5g}" if isinstance(v, float) else str(v) for v in t) + "\n")
    with open(out / "CALIBRATION.tsv", "w") as fh:
        fh.write("escenario\tX\tf_ND\tf_N\tf_D\tcociente_fN_fD\n")
        for t in calib:
            fh.write("\t".join(f"{v:.4g}" if isinstance(v, float) else str(v) for v in t) + "\n")
    print(f"{label}: {len(rows)} réplicas; sitios por rama ND/N/D: "
          f"{int(br['ND'].sum()):,}/{int(br['N'].sum()):,}/{int(br['D'].sum()):,}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]) if len(sys.argv) > 4 else 100)
