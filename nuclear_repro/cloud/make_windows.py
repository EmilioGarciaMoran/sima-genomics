#!/usr/bin/env python3
"""Define (o amplía) las ventanas de coordenadas para las distancias por pares.

Cada ventana es una unidad de trabajo independiente: se procesa en la nube por separado,
sus resultados no cambian al añadir otras, y en el jackknife cada ventana es un bloque.
Añadir ventanas solo añade líneas a pairwise/windows.bed; las existentes no se tocan.

    # una ventana por cromosoma, de 100.000 a 200.000
    python3 cloud/make_windows.py --fai hg19.fa.fai --start 100000 --size 100000
    # más adelante: ampliar con una ventana de 100 kb cada 10 Mb (se acumulan con las anteriores)
    python3 cloud/make_windows.py --fai hg19.fa.fai --start 5000000 --size 100000 --every 10000000
    # con --ref se descartan las ventanas que en la referencia son casi todo N
    # (p. ej. los primeros Mb de los cromosomas acrocéntricos 13, 14, 15, 21 y 22)

Formato: BED (0-based) chrom inicio fin id, con id = chrom_inicio1based_fin, que nombra
el directorio pairwise/windows/<id>/ con el panel y los vectores de esa ventana.
"""
import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def window_id(chrom, start0, end):
    return f"{chrom}_{start0 + 1:09d}_{end:09d}"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fai", required=True, type=Path, help="índice .fai de la referencia (longitudes)")
    ap.add_argument("--chroms", nargs="*", default=[f"chr{i}" for i in range(1, 23)])
    ap.add_argument("--start", type=int, required=True, help="inicio de la primera ventana (0-based)")
    ap.add_argument("--size", type=int, required=True)
    ap.add_argument("--every", type=int, default=0, help="repetir una ventana cada N pb (0 = solo una)")
    ap.add_argument("--ref", type=Path, help="FASTA para descartar ventanas con más de --max-n de N")
    ap.add_argument("--max-n", type=float, default=0.5)
    ap.add_argument("--out", type=Path, default=ROOT / "pairwise" / "windows.bed")
    a = ap.parse_args()

    lengths = {}
    for line in a.fai.read_text().splitlines():
        f = line.split("\t")
        lengths[f[0]] = int(f[1])
    fa = None
    if a.ref:
        import pysam
        fa = pysam.FastaFile(str(a.ref))

    existing = {}
    if a.out.exists():
        for line in a.out.read_text().splitlines():
            if line.strip() and not line.startswith("#"):
                c, s, e, i = line.split("\t")[:4]
                existing[i] = (c, int(s), int(e))
    added, skipped = [], []
    for c in a.chroms:
        if c not in lengths:
            print(f"  aviso: {c} no está en {a.fai}", file=sys.stderr)
            continue
        s = a.start
        while s + a.size <= lengths[c]:
            wid = window_id(c, s, s + a.size)
            if wid not in existing:
                if fa is not None:
                    seq = fa.fetch(c, s, s + a.size).upper()
                    if seq.count("N") / len(seq) > a.max_n:
                        skipped.append(wid)
                        s = s + a.every if a.every else lengths[c]
                        continue
                existing[wid] = (c, s, s + a.size)
                added.append(wid)
            s = s + a.every if a.every else lengths[c]

    # solapamientos entre ventanas: duplicarían sitios en la matriz
    by_chrom = {}
    for wid, (c, s, e) in existing.items():
        by_chrom.setdefault(c, []).append((s, e, wid))
    for c, ivs in by_chrom.items():
        ivs.sort()
        for (s1, e1, w1), (s2, e2, w2) in zip(ivs, ivs[1:]):
            if s2 < e1:
                sys.exit(f"ERROR: {w1} y {w2} se solapan; elige otro --start/--every")

    order = {c: i for i, c in enumerate(lengths)}
    rows = sorted(existing.items(), key=lambda kv: (order.get(kv[1][0], 1e9), kv[1][1]))
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with open(a.out, "w") as fh:
        fh.write("# ventanas de las distancias por pares (BED 0-based): chrom inicio fin id\n")
        for wid, (c, s, e) in rows:
            fh.write(f"{c}\t{s}\t{e}\t{wid}\n")
    total = sum(e - s for c, s, e in existing.values())
    print(f"{len(added)} ventanas nuevas, {len(existing)} en total ({total / 1e6:.2f} Mb) -> {a.out}")
    if skipped:
        print(f"Descartadas por N en la referencia ({len(skipped)}): " + ", ".join(skipped))


if __name__ == "__main__":
    main()
