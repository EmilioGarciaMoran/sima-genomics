#!/usr/bin/env python3
"""Filtra y elimina duplicados de un SAM de lecturas sencillas de ADN antiguo.

Lee SAM por la entrada estándar y escribe SAM por la salida estándar. Las estadísticas
van a la salida de error, para no mezclarlas con los datos.

Criterios, fijados en ../PREREGISTRO_D_statistics.md, apartado 10 ter:
  - se descartan las no mapeadas, las secundarias y las suplementarias
  - calidad de mapeo >= 30 (literal de Meyer 2016: 'a map quality score of 30
    or greater'; el preregistro decia 25, tomado por error del mq25 de los VCF)
  - longitud de la lectura >= 35 pb
  - duplicados: se conserva UNA lectura por (cromosoma, inicio, fin, hebra)

Sobre los duplicados: las librerías de ADN antiguo están muy amplificadas, así que el mismo
fragmento aparece muchas veces. La clave incluye el FIN además del inicio porque dos
fragmentos distintos pueden empezar en la misma posición y tener longitudes diferentes;
agrupar solo por inicio borraría lecturas legítimas. El fin se calcula sobre la referencia a
partir del CIGAR, no sumando la longitud de la secuencia, porque las inserciones y
deleciones hacen que ambas cosas no coincidan.

uso:
    samtools view -h entrada.bam | python3 filter_dedup.py > salida.sam
    python3 filter_dedup.py --min-mapq 37 < entrada.sam > salida.sam
"""
from __future__ import annotations

import argparse
import re
import sys

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
CONSUME_REF = frozenset("MDN=X")


def ref_len(cigar):
    """Longitud que la lectura ocupa sobre la referencia, según el CIGAR."""
    if cigar == "*":
        return 0
    return sum(int(n) for n, op in CIGAR_RE.findall(cigar) if op in CONSUME_REF)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-mapq", type=int, default=30)
    ap.add_argument("--min-len", type=int, default=35)
    ap.add_argument("--no-dedup", action="store_true",
                    help="solo para diagnóstico: no elimina duplicados")
    ap.add_argument("--etiqueta", default="")
    a = ap.parse_args()

    vistos = set()
    n = dict(total=0, sin_mapear=0, secundaria=0, mapq=0, corta=0, duplicada=0, conservada=0)
    out = sys.stdout

    for linea in sys.stdin:
        if linea.startswith("@"):
            out.write(linea)
            continue
        n["total"] += 1
        c = linea.rstrip("\n").split("\t")
        if len(c) < 11:
            continue
        flag = int(c[1])
        if flag & 0x4:
            n["sin_mapear"] += 1
            continue
        if flag & 0x100 or flag & 0x800:
            n["secundaria"] += 1
            continue
        if int(c[4]) < a.min_mapq:
            n["mapq"] += 1
            continue
        seq = c[9]
        if len(seq) < a.min_len:
            n["corta"] += 1
            continue
        if not a.no_dedup:
            ini = int(c[3])
            clave = (c[2], ini, ini + ref_len(c[5]), bool(flag & 0x10))
            if clave in vistos:
                n["duplicada"] += 1
                continue
            vistos.add(clave)
        n["conservada"] += 1
        out.write(linea)

    e = sys.stderr
    etq = f" [{a.etiqueta}]" if a.etiqueta else ""
    print(f"filter_dedup{etq}: mapq>={a.min_mapq}, longitud>={a.min_len}"
          f"{', SIN eliminar duplicados' if a.no_dedup else ''}", file=e)
    print(f"  lecturas leídas:      {n['total']:>10,}", file=e)
    print(f"  sin mapear:           {n['sin_mapear']:>10,}", file=e)
    print(f"  secundarias:          {n['secundaria']:>10,}", file=e)
    print(f"  calidad insuficiente: {n['mapq']:>10,}", file=e)
    print(f"  demasiado cortas:     {n['corta']:>10,}", file=e)
    print(f"  duplicadas:           {n['duplicada']:>10,}", file=e)
    print(f"  CONSERVADAS:          {n['conservada']:>10,}", file=e)
    if n["total"]:
        print(f"  conservadas / leídas: {100.0*n['conservada']/n['total']:>9.1f}%", file=e)


if __name__ == "__main__":
    main()
