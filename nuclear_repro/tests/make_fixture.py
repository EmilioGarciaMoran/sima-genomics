#!/usr/bin/env python3
"""Genera un caso sintético con resultados conocidos para probar todo el flujo sin datos reales.

Genoma de 2 cromosomas (chr1 de 123.456 pb -> el ensamblaje se detecta como 'fixture'),
dos muestras (SH_A, SH_B) con lecturas de 40 pb, y verdades plantadas:

  Catálogo                                     Esperado
  P1 chr1:50001  REF=der, 10 lecturas der      SH_A: DER
  P2 chr1:80001  4+2 lecturas ancestrales      SH_A, SH_B, POOL: ANC
  P3 chr1:110001 lecturas a 1 kb, ninguna encima  NA, reads_window > 0
  P4 chr2:30001  nada en ± 2 kb                NA, reads_window = 0
  P5 chr2:60001  C>T, 3 lecturas con T en la posición 0 (daño)
                                               DER (transición) pero NA sin extremos
  CDS
  GENE_HOM chr1:[90000,90300)  A->C homocigoto en 90101, 8 lecturas  -> CC
  GENE_HET chr1:[100000,100200) G/T heterocigoto en 100051, 6+6     -> GT
  GENE_EMPTY chr2:[100000,100150) sin lecturas                      -> 0 pb cubiertos
  Lecturas fuera de toda región (chr1:118000) -> el recorte las tiene que eliminar.

    python3 tests/make_fixture.py OUTDIR   ->  OUTDIR/full/{ref.fa,SH_A.bam,SH_B.bam}, OUTDIR/catalog/
"""
import random
import sys
from pathlib import Path

import pysam

LEN = {"chr1": 123456, "chr2": 120000}
RL = 40
rng = random.Random(12345)


def make_ref():
    ref = {c: [rng.choice("ACGT") for _ in range(n)] for c, n in LEN.items()}
    fixed = {("chr1", 50001): "C", ("chr1", 80001): "G", ("chr1", 110001): "A", ("chr2", 30001): "T",
             ("chr2", 60001): "C", ("chr1", 90101): "A", ("chr1", 100051): "G"}
    for (c, p), b in fixed.items():
        ref[c][p - 1] = b
    return {c: "".join(s) for c, s in ref.items()}


def read(header, ref, name, chrom, start0, reverse=False, subs=None):
    seq = list(ref[chrom][start0:start0 + RL])
    for pos1, b in (subs or {}).items():
        if start0 <= pos1 - 1 < start0 + RL:
            seq[pos1 - 1 - start0] = b
    a = pysam.AlignedSegment(header)
    a.query_name, a.reference_name, a.reference_start = name, chrom, start0
    a.flag, a.mapping_quality, a.cigarstring = (16 if reverse else 0), 37, f"{RL}M"
    a.query_sequence = "".join(seq)
    a.query_qualities = pysam.qualitystring_to_array("?" * RL)  # Q30
    return a


def covering(pos1, n, offset=7):
    """n posiciones de inicio (0-based) que cubren pos1 sin caer en los extremos de la lectura."""
    return [pos1 - 1 - offset - (i % 20) for i in range(n)]


def main(out: Path):
    ref = make_ref()
    (out / "full").mkdir(parents=True, exist_ok=True)
    (out / "catalog").mkdir(parents=True, exist_ok=True)
    with open(out / "full" / "ref.fa", "w") as fh:
        for c, s in ref.items():
            fh.write(f">{c}\n" + "\n".join(s[i:i + 60] for i in range(0, len(s), 60)) + "\n")
    pysam.faidx(str(out / "full" / "ref.fa"))
    header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "unsorted"},
        "SQ": [{"SN": c, "LN": n} for c, n in LEN.items()]})

    plan = {"SH_A": [], "SH_B": []}
    k = 0

    def add(sample, chrom, starts, subs=None):
        nonlocal k
        for s in starts:
            plan[sample].append(read(header, ref, f"r{k}", chrom, s, reverse=k % 2 == 1, subs=subs))
            k += 1

    add("SH_A", "chr1", covering(50001, 10))
    add("SH_A", "chr1", covering(80001, 4))
    add("SH_B", "chr1", covering(80001, 2))
    add("SH_A", "chr1", [111000, 111050])                                # cerca de P3, sin cubrirla
    for _ in range(3):                                                   # daño: T en la posición 0
        a = read(header, ref, f"r{k}", "chr2", 60000, subs={60001: "T"}); k += 1
        plan["SH_A"].append(a)
    add("SH_A", "chr1", covering(90101, 8), subs={90101: "C"})           # homocigoto alternativo
    add("SH_A", "chr1", covering(100051, 6), subs={100051: "T"})         # heterocigoto
    add("SH_A", "chr1", covering(100051, 6, offset=9))
    add("SH_B", "chr1", [118000, 118010, 118020])                        # fuera de toda región

    for sample, reads in plan.items():
        tmp = out / "full" / f"{sample}.unsorted.bam"
        with pysam.AlignmentFile(str(tmp), "wb", header=header) as fh:
            for r in reads:
                fh.write(r)
        pysam.sort("-o", str(out / "full" / f"{sample}.bam"), str(tmp))
        tmp.unlink()
        pysam.index(str(out / "full" / f"{sample}.bam"))

    (out / "catalog" / "positions.tsv").write_text(
        "## catálogo sintético\n"
        "id\tchrom\tpos\tref\tanc\tder\tgene\tpanel\tsource\n"
        "P1\tchr1\t50001\tC\tA\tC\tGENE1\ttest\tfixture\n"
        "P2\tchr1\t80001\tG\tG\tT\tGENE2\ttest\tfixture\n"
        "P3\tchr1\t110001\tA\tA\tG\tGENE3\ttest\tfixture\n"
        "P4\tchr2\t30001\tT\tT\tG\tGENE4\ttest\tfixture\n"
        "P5\tchr2\t60001\tC\tC\tT\tGENE5\ttest\tfixture\n")
    (out / "catalog" / "cds_regions.bed").write_text(
        "# CDS sintéticas\n"
        "chr1\t90000\t90300\tGENE_HOM\n"
        "chr1\t100000\t100200\tGENE_HET\n"
        "chr2\t100000\t100150\tGENE_EMPTY\n")
    print(f"fixture en {out}")


if __name__ == "__main__":
    main(Path(sys.argv[1]))
