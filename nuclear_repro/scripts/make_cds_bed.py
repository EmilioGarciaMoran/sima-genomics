#!/usr/bin/env python3
"""Genera catalog/cds_regions.bed a partir de un GTF de Ensembl y una lista de genes.

Toma las filas 'CDS' de los genes pedidos, las pasa a BED 0-based y fusiona los
solapamientos entre isoformas (equivale a bedtools merge por gen).

    python3 scripts/make_cds_bed.py --gtf Homo_sapiens.GRCh37.87.gtf.gz \
        --genes FOXP2 CNTNAP2 SCN9A SLC24A5 MC1R SLC45A2 TKTL1 NOVA1 ASPM TLR1 OAS1 ADSL BAZ1B ... \
        [--ucsc-chroms] > catalog/cds_regions.bed
"""
import argparse
import gzip
import re
import sys
from collections import defaultdict


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--genes", required=True, nargs="+")
    ap.add_argument("--ucsc-chroms", action="store_true", help="añadir prefijo 'chr' (Ensembl usa '7', UCSC 'chr7')")
    a = ap.parse_args()

    wanted = set(a.genes)
    raw = defaultdict(list)
    opener = gzip.open if a.gtf.endswith(".gz") else open
    with opener(a.gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[2] != "CDS":
                continue
            m = re.search(r'gene_name "([^"]+)"', f[8])
            if not m or m.group(1) not in wanted:
                continue
            chrom = ("chr" + f[0]) if a.ucsc_chroms and not f[0].startswith("chr") else f[0]
            raw[m.group(1)].append((chrom, int(f[3]) - 1, int(f[4])))

    n_raw = n_merged = bp = 0
    out = []
    for gene, ivs in raw.items():
        n_raw += len(ivs)
        merged = []
        for c, s, e in sorted(ivs):
            if merged and merged[-1][0] == c and s <= merged[-1][2]:
                merged[-1][2] = max(merged[-1][2], e)
            else:
                merged.append([c, s, e])
        for c, s, e in merged:
            out.append((c, s, e, gene))
            bp += e - s
        n_merged += len(merged)
    for row in sorted(out):
        print("\t".join(map(str, row)))
    absent = sorted(wanted - set(raw))
    print(f"{n_raw} CDS -> {n_merged} regiones fusionadas, {bp:,} pb"
          + (f"; genes no encontrados: {', '.join(absent)}" if absent else ""), file=sys.stderr)


if __name__ == "__main__":
    main()
