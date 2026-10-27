#!/usr/bin/env python3
"""Pieza 4: recuento por rama y D-statistics con los alelos de Sima.

Aplica lo fijado en ../PREREGISTRO_D_statistics.md. Produce tres cosas:

1. REPRODUCCIÓN DE MEYER 2016. Fracción de alelos derivados que Sima comparte con cada
   rama, definidas con Denisova 5 (Altai) y Denisova 3, igual que el artículo:
       rama ND  = ambos derivados          Meyer: 87% (AT-5431), 68% (AT-5482)
       rama N   = solo Denisova 5          Meyer: 43%, 39%
       rama D   = solo Denisova 3          Meyer:  9%,  7%
   Es el control de validación del pipeline: si no se parecen, algo está mal.

2. D-STATISTIC PREREGISTRADA. D(Denisova 3, Vindija 33.19; Sima, ancestral).
       D > 0  -> Sima comparte más con Vindija (linaje neandertal)  = H1
       D ≈ 0  -> Sima equidistante de ambos linajes                 = H2 (Homo atapuercae)
   Error típico por jackknife de bloques ponderado (Busing et al. 1999), bloques de 5 Mb.

3. PRUEBA DE EQUIVALENCIA. Con margen δ = 0,10: H2 solo se sostiene si el intervalo de
   confianza del 90% cabe entero dentro de −δ a +δ.

uso:
    python3 scripts/analyze_sima.py --sitios ../branch_counts/results/genoma_inf.npz \\
            --alelos results/alelos --out results/RESULTADO.txt
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from scipy.stats import binomtest, norm

# branch_counts.py vive en ../../branch_counts/scripts/ : scripts -> sima_pipeline ->
# nuclear_repro -> branch_counts/scripts. Se prueban varias rutas por si el script se
# invoca desde otro directorio.
_aqui = Path(__file__).resolve().parent
for _c in (_aqui.parent.parent / "branch_counts" / "scripts",
           _aqui.parent / "branch_counts" / "scripts",
           _aqui):
    if (_c / "branch_counts.py").exists():
        sys.path.insert(0, str(_c))
        break
else:
    sys.exit("ERROR: no encuentro branch_counts.py; se esperaba en "
             f"{_aqui.parent.parent / 'branch_counts' / 'scripts'}")
from branch_counts import branch_counts  # noqa: E402

DELTA = 0.10      # margen de equivalencia, preregistrado
ZCRIT = 3.0       # umbral de significación, preregistrado
BLOCK = 5_000_000


def jackknife_D(abba_b, baba_b):
    """D y su error típico por jackknife de bloques ponderado (Busing et al. 1999)."""
    n_b = abba_b + baba_b
    A, B, N = abba_b.sum(), baba_b.sum(), n_b.sum()
    if N == 0:
        return np.nan, np.nan, 0
    D = (A - B) / N
    keep = n_b > 0
    if keep.sum() < 2:
        return D, np.nan, int(N)
    ab, bb, nb = abba_b[keep], baba_b[keep], n_b[keep]
    g = len(nb)
    with np.errstate(invalid="ignore", divide="ignore"):
        Dj = np.where(N - nb > 0, ((A - ab) - (B - bb)) / (N - nb), np.nan)
    h = N / nb
    tau = h * D - (h - 1) * Dj
    ok = np.isfinite(tau)
    if ok.sum() < 2:
        return D, np.nan, int(N)
    tau, h = tau[ok], h[ok]
    var = np.sum((tau - tau.mean()) ** 2 / (h - 1)) / g
    return D, (float(np.sqrt(var)) if var > 0 else np.nan), int(N)


def calcula_D(sima, d3, vi, mbuti, blocks):
    """D(D3, Vindija; Sima, ancestral) a partir de los códigos 0/1/-1."""
    base = (d3 >= 0) & (vi >= 0) & (mbuti == 0) & (sima == 1)
    abba = base & (vi == 1) & (d3 == 0)      # Sima con Vindija
    baba = base & (d3 == 1) & (vi == 0)      # Sima con Denisova
    nb = blocks.max() + 1
    return jackknife_D(np.bincount(blocks[abba], minlength=nb).astype(float),
                       np.bincount(blocks[baba], minlength=nb).astype(float)), \
           int(abba.sum()), int(baba.sum())


def linea(f, txt=""):
    print(txt)
    f.write(txt + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sitios", required=True)
    ap.add_argument("--alelos", required=True)
    ap.add_argument("--out", default="results/RESULTADO.txt")
    a = ap.parse_args()

    z = np.load(a.sitios, allow_pickle=False)
    pos, chrom = z["pos"], z["chrom"]
    D5, VI, D3, CH, MB = z["D5"], z["Vi3319"], z["D3"], z["Chag8"], z["Mbuti"]
    blk = chrom.astype(np.int64) * 10_000 + (pos // BLOCK)
    _, blocks = np.unique(blk, return_inverse=True)

    alel = Path(a.alelos)
    def carga(nombre):
        f = alel / f"{nombre}.npz"
        return np.load(f, allow_pickle=False) if f.exists() else None

    out = Path(a.out); out.parent.mkdir(parents=True, exist_ok=True)
    fh = open(out, "w")
    L = lambda t="": linea(fh, t)  # noqa: E731

    L("=" * 70)
    L("RESULTADO DEL ANÁLISIS NUCLEAR DE SIMA DE LOS HUESOS")
    L("=" * 70)
    L(f"Sitios de referencia: {len(pos):,}   bloques de jackknife: {blocks.max()+1}")
    L("Parámetros y criterios: ../PREREGISTRO_D_statistics.md, fijados antes de los datos.")
    L()

    # ---------- 1. reproducción de Meyer ----------
    L("-" * 70)
    L("1. REPRODUCCIÓN DE MEYER 2016 (control de validación)")
    L("-" * 70)
    L("   Ramas definidas con Denisova 5 (Altai) y Denisova 3, como en el artículo.")
    L()
    L(f"   {'espécimen':20s} {'n_ND':>7s} {'f_ND':>7s} {'Meyer':>6s} | "
      f"{'n_N':>6s} {'f_N':>6s} {'Meyer':>6s} | {'n_D':>6s} {'f_D':>6s} {'Meyer':>6s}")
    MEYER = {"AT-5431_femur": (0.87, 0.43, 0.09), "AT-5482_incisivo": (0.68, 0.39, 0.07)}
    for esp in ("AT-5431_femur", "AT-5482_incisivo", "femurXIII", "AT-5444_molar",
                "AT-6672_escapula"):
        d = carga(esp)
        if d is None:
            continue
        r = branch_counts(d["codigo"], D5, D3, MB, blocks=blocks)
        m = MEYER.get(esp, (np.nan,) * 3)
        def c(k, i):
            v = r[k]
            f = f"{v['f']:.2f}" if v["n"] else "  -"
            mm = f"{m[i]:.2f}" if np.isfinite(m[i]) else "  -"
            return f"{v['n']:>6,} {f:>6s} {mm:>6s}"
        L(f"   {esp:20s} {c('AB',0)} | {c('A',1)} | {c('B',2)}")
    L()
    L("   f_ND alto y f_D bajo indica pertenencia al linaje neandertal.")
    L("   Nota: ENA contiene ~60% de las lecturas de Meyer, de modo que las n son menores;")
    L("   las fracciones sí son comparables.")
    L()

    # ---------- 2. D-statistic ----------
    L("-" * 70)
    L("2. D-STATISTIC PREREGISTRADA: D(Denisova 3, Vindija 33.19; Sima, ancestral)")
    L("-" * 70)
    L("   D > 0: Sima comparte más con Vindija (H1, linaje neandertal)")
    L("   D ≈ 0: Sima equidistante (H2, linaje aparte / Homo atapuercae)")
    L()
    L(f"   {'espécimen':22s} {'ABBA':>6s} {'BABA':>6s} {'n':>6s} {'D':>8s} {'EE':>7s} "
      f"{'Z':>7s} {'IC90':>16s}  decisión")
    resultados = {}
    combinado = {}
    for esp in ("AT-5482_incisivo", "AT-5431_femur", "femurXIII", "AT-5444_molar",
                "AT-6672_escapula"):
        d = carga(esp)
        if d is None:
            continue
        (D, ee, n), ab, ba = calcula_D(d["codigo"], D3, VI, MB, blocks)
        resultados[esp] = (D, ee, n, ab, ba)
        if esp in ("AT-5482_incisivo", "AT-5431_femur"):
            combinado[esp] = d["codigo"]
        if not n:
            L(f"   {esp:22s} {'sin sitios informativos':>50s}")
            continue
        Zv = D / ee if ee and np.isfinite(ee) else np.nan
        lo, hi = (D - 1.645 * ee, D + 1.645 * ee) if np.isfinite(ee) else (np.nan, np.nan)
        dec = ("H2 DESCARTADA" if np.isfinite(Zv) and abs(Zv) > ZCRIT else
               ("H2 sostenible" if np.isfinite(lo) and -DELTA <= lo and hi <= DELTA
                else "no concluyente"))
        L(f"   {esp:22s} {ab:>6,} {ba:>6,} {n:>6,} {D:>8.3f} {ee:>7.3f} {Zv:>7.1f} "
          f"[{lo:>6.2f},{hi:>6.2f}]  {dec}")

    # combinado de los dos principales
    if len(combinado) == 2:
        c1, c2 = combinado["AT-5482_incisivo"], combinado["AT-5431_femur"]
        # un sitio cubierto por ambos: se toma el del incisivo, que tiene más datos
        comb = np.where(c1 >= 0, c1, c2)
        (D, ee, n), ab, ba = calcula_D(comb, D3, VI, MB, blocks)
        resultados["COMBINADO"] = (D, ee, n, ab, ba)
        if n:
            Zv = D / ee if ee and np.isfinite(ee) else np.nan
            lo, hi = (D - 1.645 * ee, D + 1.645 * ee) if np.isfinite(ee) else (np.nan,) * 2
            dec = ("H2 DESCARTADA" if np.isfinite(Zv) and abs(Zv) > ZCRIT else
                   ("H2 sostenible" if np.isfinite(lo) and -DELTA <= lo and hi <= DELTA
                    else "no concluyente"))
            L(f"   {'>> COMBINADO (principal)':22s} {ab:>6,} {ba:>6,} {n:>6,} {D:>8.3f} "
              f"{ee:>7.3f} {Zv:>7.1f} [{lo:>6.2f},{hi:>6.2f}]  {dec}")
            p = binomtest(ab, n, 0.5).pvalue if n else np.nan
            L()
            L(f"   Prueba binomial exacta del combinado: p = {p:.3g}")
            L(f"   (se incluye porque con pocos sitios la aproximación normal no es fiable)")
    L()

    # ---------- 3. reglas de interpretación ----------
    L("-" * 70)
    L("3. REGLAS DE INTERPRETACIÓN PREREGISTRADAS (apartado 9 del preregistro)")
    L("-" * 70)
    if "COMBINADO" in resultados:
        D = resultados["COMBINADO"][0]
        if D >= 0.50:
            v = ("D >= 0,50: CONFIRMA MEYER 2016. Sima dentro del linaje neandertal.\n"
                 "   H2 (Homo atapuercae) descartada por la vía genética.")
        elif D >= 0.13:
            v = ("D entre 0,13 y 0,50: significativa pero MUY POR DEBAJO de lo predicho.\n"
                 "   NO INTERPRETAR BIOLÓGICAMENTE. Investigar antes contaminación, sesgo\n"
                 "   hacia la referencia y la traducción de los porcentajes de Meyer.")
        elif D > -0.13:
            v = ("D < 0,13: incompatible con la predicción derivada de Meyer.\n"
                 "   Revisar razonamiento y datos antes de reclamar nada.")
        else:
            v = ("D NEGATIVA: Sima compartiría más con Denisova. No lo predice ninguna\n"
                 "   hipótesis. Error técnico casi seguro.")
        L(f"   {v}")
    L()

    # ---------- 4. sensibilidad ----------
    L("-" * 70)
    L("4. ANÁLISIS DE SENSIBILIDAD PREREGISTRADOS")
    L("-" * 70)
    L(f"   {'variante':30s} {'n':>6s} {'D':>8s} {'EE':>7s} {'Z':>7s}")
    for esp in ("AT-5482_incisivo", "AT-5431_femur", "femurXIII"):
        for suf, etq in (("", "principal"), (".nofilt", "sin filtro de daño"),
                         (".ss", "criterio cadena sencilla"), (".ds", "criterio doble cadena")):
            d = carga(esp + suf)
            if d is None:
                continue
            (D, ee, n), _, _ = calcula_D(d["codigo"], D3, VI, MB, blocks)
            if not n:
                continue
            Zv = D / ee if ee and np.isfinite(ee) else np.nan
            L(f"   {esp.split('_')[0]+' '+etq:30s} {n:>6,} {D:>8.3f} {ee:>7.3f} {Zv:>7.1f}")
    L()
    L("   Un resultado que solo aparezca en la variante principal y se desvanezca en las")
    L("   tres sensibilidades NO se reporta como hallazgo (preregistro, apartado 10 ter).")
    L()

    # ---------- 5. dispersión del muestreo ----------
    L("-" * 70)
    L("5. DISPERSIÓN POR EL MUESTREO PSEUDOHAPLOIDE")
    L("-" * 70)
    for esp in ("AT-5482_incisivo", "AT-5431_femur"):
        d = carga(esp)
        if d is None or "codigos" not in d.files:
            continue
        Ds = []
        for r in range(d["codigos"].shape[0]):
            (D, _, n), _, _ = calcula_D(d["codigos"][r], D3, VI, MB, blocks)
            if n:
                Ds.append(D)
        if len(Ds) > 1:
            ee = resultados[esp][1]
            L(f"   {esp:22s} D entre {min(Ds):.3f} y {max(Ds):.3f}, "
              f"desviación {np.std(Ds):.4f}")
            if np.isfinite(ee) and ee > 0:
                L(f"   {'':22s} frente al EE del jackknife {ee:.3f}: "
                  f"{'despreciable' if np.std(Ds) < ee/5 else 'NO despreciable'}")
    L()
    L("=" * 70)
    fh.close()
    # La ruta NO se escribe dentro del fichero: si lo hiciera, el propio nombre de salida
    # cambiaría el contenido y la comparación byte a byte dejaría de ser posible.
    print(f"Escrito en {out}")


if __name__ == "__main__":
    main()
