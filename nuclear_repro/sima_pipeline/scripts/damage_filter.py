#!/usr/bin/env python3
"""Pieza 2: filtro de daño y enmascarado terminal, al modo de Meyer et al. 2016.

Lee SAM por la entrada estándar y escribe SAM por la salida estándar. Las estadísticas van
a la salida de error.

QUÉ HACE
  1. Conserva solo los fragmentos que muestran desaminación en un extremo (daño terminal).
     Es lo que enriquece en ADN endógeno frente a la contaminación moderna, que no lo tiene.
  2. En los conservados, pone a N las bases de los 3 primeros y 3 últimos lugares alineados
     que podrían proceder de daño, para que el propio daño no genere alelos falsos.

LIBRERÍAS DE CADENA SENCILLA (--library ss, por defecto)
  Meyer 2016 usó el protocolo de cadena sencilla de Gansauge y Meyer, que es lo que permitió
  recuperar ADN tan degradado. En esas librerías la desaminación aparece como C->T en AMBOS
  extremos de la cadena original. Al alinear:
     - lectura en hebra directa: se ve C->T en los dos extremos
     - lectura en hebra inversa: el SAM guarda el complemento inverso, así que se ve G->A
  Por eso el criterio y el enmascarado dependen de la hebra, que se lee del FLAG.
  Esto explica la frase del artículo "T enmascaradas en las 3 primeras y 3 últimas
  posiciones": está escrita en la orientación de la cadena original, no del SAM.

LIBRERÍAS DE DOBLE CADENA (--library ds)
  Convención alternativa, para comparar: C->T en el extremo 5' y G->A en el 3', en
  orientación de referencia, independientemente de la hebra.

POSICIONES
  Se cuentan sobre la porción ALINEADA, descartando el recorte blando: las bases recortadas
  no se alinean, así que no contribuyen a ninguna llamada de alelo.

uso:
    samtools view -h entrada.bam | python3 damage_filter.py --ref ref.fasta > salida.sam
    ... --library ds          convención de doble cadena
    ... --mask 0              no enmascarar (solo para diagnóstico)
    ... --no-filter           no filtrar por daño (análisis de sensibilidad preregistrado)
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from fasta import Fasta  # noqa: E402

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
CONSUME_REF = frozenset("MDN=X")
CONSUME_READ = frozenset("MIS=X")


def parse_cigar(cigar):
    return [(int(n), op) for n, op in CIGAR_RE.findall(cigar)] if cigar != "*" else []


def aligned_pairs(cigar, pos0):
    """Devuelve [(indice_en_SEQ, posicion_en_referencia)] de las bases alineadas.

    pos0 es 0-based. Se omiten inserciones y recorte blando: no tienen base de referencia.
    """
    pares = []
    r = pos0      # posición en la referencia
    q = 0         # índice en SEQ
    for n, op in parse_cigar(cigar):
        if op in ("M", "=", "X"):
            for i in range(n):
                pares.append((q + i, r + i))
            q += n; r += n
        elif op in ("I", "S"):
            q += n
        elif op in ("D", "N"):
            r += n
        # H y P no consumen SEQ
    return pares


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True)
    ap.add_argument("--library", choices=["ss", "ds"], default="ss")
    ap.add_argument("--damage-window", type=int, default=1,
                    help="en cuántas posiciones terminales se busca el daño. El artículo "
                         "dice 'C->T terminal' sin precisar; 1 es la lectura estricta y 3 la "
                         "convención habitual del MPI.")
    ap.add_argument("--mask", type=int, default=3,
                    help="cuántas posiciones enmascarar en cada extremo (0 = ninguna)")
    ap.add_argument("--no-filter", action="store_true",
                    help="no exige daño terminal (análisis de sensibilidad preregistrado)")
    ap.add_argument("--etiqueta", default="")
    a = ap.parse_args()

    ref = Fasta(a.ref)
    out = sys.stdout
    n = dict(total=0, sin_alinear=0, sin_dano=0, conservada=0, bases=0, enmascaradas=0,
             dano_5=0, dano_3=0, dano_ambos=0)

    for linea in sys.stdin:
        if linea.startswith("@"):
            out.write(linea)
            continue
        n["total"] += 1
        c = linea.rstrip("\n").split("\t")
        if len(c) < 11:
            continue
        flag = int(c[1])
        chrom, pos0, cigar, seq = c[2], int(c[3]) - 1, c[5], c[9]
        pares = aligned_pairs(cigar, pos0)
        if not pares:
            n["sin_alinear"] += 1
            continue

        rini, rfin = pares[0][1], pares[-1][1]
        refseq = ref.fetch(chrom, rini, rfin + 1)
        if len(refseq) != rfin - rini + 1:
            n["sin_alinear"] += 1
            continue

        def base_ref(rpos):
            return refseq[rpos - rini]

        inversa = bool(flag & 0x10)
        # Qué sustitución delata daño en cada extremo, según la librería y la hebra.
        if a.library == "ss":
            # cadena sencilla: C->T en ambos extremos de la cadena original
            izq = der = ("G", "A") if inversa else ("C", "T")
        else:
            # doble cadena: C->T en el 5' y G->A en el 3', en orientación de referencia
            izq, der = ("C", "T"), ("G", "A")

        w = max(1, a.damage_window)
        d5 = any(base_ref(r) == izq[0] and seq[q].upper() == izq[1] for q, r in pares[:w])
        d3 = any(base_ref(r) == der[0] and seq[q].upper() == der[1] for q, r in pares[-w:])

        if not (d5 or d3):
            if not a.no_filter:
                n["sin_dano"] += 1
                continue
        else:
            n["dano_5"] += d5 and not d3
            n["dano_3"] += d3 and not d5
            n["dano_ambos"] += d5 and d3

        # Enmascarado: las bases que podrían venir de daño, en los m primeros y m últimos
        # lugares alineados, se ponen a N.
        if a.mask > 0:
            s = list(seq)
            m = a.mask
            sospechosa_izq = izq[1]     # la base RESULTANTE del daño (T o A)
            sospechosa_der = der[1]
            for q, _ in pares[:m]:
                if s[q].upper() == sospechosa_izq:
                    s[q] = "N"; n["enmascaradas"] += 1
            for q, _ in pares[-m:]:
                if s[q].upper() == sospechosa_der:
                    s[q] = "N"; n["enmascaradas"] += 1
            c[9] = "".join(s)

        n["conservada"] += 1
        n["bases"] += len(pares)
        out.write("\t".join(c) + "\n")

    e = sys.stderr
    etq = f" [{a.etiqueta}]" if a.etiqueta else ""
    print(f"damage_filter{etq}: librería {a.library}, ventana de daño {a.damage_window}, "
          f"enmascarado {a.mask}+{a.mask}"
          f"{', SIN filtrar por daño' if a.no_filter else ''}", file=e)
    print(f"  lecturas leídas:    {n['total']:>10,}", file=e)
    print(f"  sin alineación útil:{n['sin_alinear']:>10,}", file=e)
    if not a.no_filter:
        print(f"  sin daño terminal:  {n['sin_dano']:>10,}", file=e)
    print(f"  CONSERVADAS:        {n['conservada']:>10,}", file=e)
    if n["total"]:
        print(f"  tasa de daño:       {100.0*n['conservada']/n['total']:>9.1f}%", file=e)
    print(f"    solo extremo izq: {n['dano_5']:>10,}", file=e)
    print(f"    solo extremo der: {n['dano_3']:>10,}", file=e)
    print(f"    ambos extremos:   {n['dano_ambos']:>10,}", file=e)
    print(f"  bases alineadas:    {n['bases']:>10,}", file=e)
    print(f"  bases enmascaradas: {n['enmascaradas']:>10,}", file=e)
    print(f"  BASES ÚTILES:       {n['bases']-n['enmascaradas']:>10,}", file=e)
    ref.close()


if __name__ == "__main__":
    main()
