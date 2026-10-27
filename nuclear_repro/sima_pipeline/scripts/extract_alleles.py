#!/usr/bin/env python3
"""Pieza 3: extrae los alelos de Sima en los sitios de referencia.

Para cada uno de los sitios informativos construidos en ../branch_counts/, mira qué bases
aporta el BAM del espécimen y codifica el resultado igual que los genomas de alta cobertura:

    -1  sin dato        0  alelo ancestral        1  alelo derivado

MUESTREO PSEUDOHAPLOIDE
  Cuando varias lecturas cubren el mismo sitio se toma UNA al azar, no un consenso. Con
  cobertura inferior a 1x el consenso no tiene sentido, y además introduciría un sesgo:
  llamar un genotipo a partir de dos lecturas favorece al alelo que coincide con la
  referencia. Es lo que exige el preregistro para el individuo de prueba, y lo que hace que
  el resultado sea comparable con la simulación de potencia.

BASES ENMASCARADAS
  Las N que dejó el filtro de daño cuentan como ausencia de dato, no como desacuerdo.

TERCER ALELO
  Si la base observada no es ni la ancestral ni la derivada del sitio, se descarta: puede
  venir de error de secuenciación o de un alelo distinto, y en ninguno de los dos casos
  informa sobre la rama.

uso:
    python3 scripts/extract_alleles.py \\
        --sitios ../branch_counts/results/genoma_inf.npz \\
        --bam ~/sima/sima_data/bam_dmg/AT-5482_incisivo.bam \\
        --etiqueta AT-5482_incisivo \\
        --out results/alelos/AT-5482_incisivo.npz
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sitios", required=True)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--etiqueta", default="")
    ap.add_argument("--seed", type=int, default=20261025)
    ap.add_argument("--repeticiones", type=int, default=20,
                    help="cuántos muestreos pseudohaploides independientes se guardan. "
                         "Fijar la semilla da reproducibilidad; repetir el muestreo mide "
                         "cuánto se mueve el resultado por el propio azar, que es una "
                         "fuente de variación distinta del jackknife. La repetición 0 es "
                         "la principal; las demás sirven para cuantificar esa dispersión.")
    ap.add_argument("--sam", action="store_true",
                    help="la entrada es un SAM de texto en vez de un BAM (para pruebas)")
    a = ap.parse_args()

    z = np.load(a.sitios, allow_pickle=False)
    pos, chrom = z["pos"], z["chrom"]
    canc, der = z["canc"], z["der"]
    nsitios = len(pos)

    # índice: (cromosoma, posición) -> fila. Las posiciones del npz son 1-based, como el SAM.
    clave = (chrom.astype(np.int64) << 32) | pos.astype(np.int64)
    indice = {int(k): i for i, k in enumerate(clave)}

    obs = defaultdict(list)     # fila -> bases observadas
    n = dict(lecturas=0, alineadas=0, en_sitio=0, enmascaradas=0, tercer_alelo=0,
             cigar_incoherente=0)

    if a.sam:
        fuente = open(a.bam)
        proc = None
    else:
        proc = subprocess.Popen(["samtools", "view", a.bam], stdout=subprocess.PIPE,
                                text=True, bufsize=1 << 20)
        fuente = proc.stdout

    for linea in fuente:
        if linea.startswith("@"):
            continue
        c = linea.rstrip("\n").split("\t")
        if len(c) < 11:
            continue
        n["lecturas"] += 1
        try:
            cr = int(c[2])
        except ValueError:
            continue                      # cromosomas no autosómicos (X, Y, MT)
        base_cr = cr << 32
        r = int(c[3])                     # 1-based
        q = 0
        seq = c[9]
        cig = CIGAR_RE.findall(c[5])
        # El CIGAR debe consumir exactamente tantas bases de la secuencia como tiene.
        # Si no, el registro está mal formado: se descarta en lugar de leer fuera de rango.
        consume = sum(int(x) for x, o in cig if o in ("M", "I", "S", "=", "X"))
        if consume != len(seq):
            n["cigar_incoherente"] += 1
            continue
        for num, op in cig:
            num = int(num)
            if op in ("M", "=", "X"):
                for i in range(num):
                    n["alineadas"] += 1
                    fila = indice.get(base_cr | (r + i))
                    if fila is None:
                        continue
                    n["en_sitio"] += 1
                    b = seq[q + i].upper()
                    if b == "N":
                        n["enmascaradas"] += 1
                        continue
                    obs[fila].append(b)
                q += num; r += num
            elif op in ("I", "S"):
                q += num
            elif op in ("D", "N"):
                r += num

    if proc is not None:
        proc.stdout.close(); proc.wait()
    else:
        fuente.close()

    # muestreo pseudohaploide y codificación
    R = max(1, a.repeticiones)
    rng = np.random.default_rng(a.seed)
    codigos = np.full((R, nsitios), -1, dtype=np.int8)
    cob = np.zeros(nsitios, dtype=np.int16)
    for fila, bases in obs.items():
        cob[fila] = min(len(bases), 32767)
        ca, de = canc[fila], der[fila]
        if len(bases) == 1:
            # un solo dato: el muestreo no interviene, el valor es el mismo en todas
            bb = ord(bases[0])
            v = 0 if bb == ca else (1 if bb == de else -1)
            if v < 0:
                n["tercer_alelo"] += 1
            codigos[:, fila] = v
        else:
            elegidas = rng.integers(len(bases), size=R)
            for r in range(R):
                bb = ord(bases[elegidas[r]])
                codigos[r, fila] = 0 if bb == ca else (1 if bb == de else -1)
            if any(ord(b) not in (ca, de) for b in bases):
                n["tercer_alelo"] += 1
    codigo = codigos[0]

    out = Path(a.out); out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out, codigo=codigo, codigos=codigos, cobertura=cob,
                        etiqueta=np.array([a.etiqueta]), seed=np.array([a.seed]),
                        repeticiones=np.array([R]), sitios=np.array([a.sitios]))

    cubiertos = int((codigo >= 0).sum())
    derivados = int((codigo == 1).sum())
    e = sys.stderr
    etq = f" [{a.etiqueta}]" if a.etiqueta else ""
    print(f"extract_alleles{etq}", file=e)
    print(f"  sitios de referencia:   {nsitios:>10,}", file=e)
    print(f"  lecturas procesadas:    {n['lecturas']:>10,}", file=e)
    print(f"  bases alineadas:        {n['alineadas']:>10,}", file=e)
    print(f"  caen en un sitio:       {n['en_sitio']:>10,}", file=e)
    print(f"    de ellas enmascaradas:{n['enmascaradas']:>10,}", file=e)
    print(f"    tercer alelo:         {n['tercer_alelo']:>10,}", file=e)
    if n["cigar_incoherente"]:
        print(f"  AVISO: {n['cigar_incoherente']:,} lectura(s) con CIGAR incoherente, "
              f"descartadas", file=e)
    print(f"  SITIOS CON DATO:        {cubiertos:>10,}", file=e)
    print(f"    derivados:            {derivados:>10,}", file=e)
    print(f"    ancestrales:          {cubiertos-derivados:>10,}", file=e)
    if cubiertos:
        rep = int((cob > 1).sum())
        print(f"  sitios con mas de una lectura: {rep:,} "
              f"({100.0*rep/cubiertos:.1f}% de los cubiertos)", file=e)
        print(f"    son los unicos donde el muestreo pseudohaploide interviene", file=e)
        if R > 1:
            d = codigos[:, codigos[0] >= 0]
            var = int((d != d[0]).any(axis=0).sum())
            print(f"  con {R} muestreos, sitios que cambian de valor: {var:,} "
                  f"({100.0*var/cubiertos:.2f}% de los cubiertos)", file=e)
    print(f"  escrito: {out}", file=e)


if __name__ == "__main__":
    main()
