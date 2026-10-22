#!/usr/bin/env python3
"""NUBE — paso 9: define el panel de sitios informativos, sin mirar a Sima.

Un sitio entra en el panel si es una transversión bialélica y algún genoma de
comparación (VCF de alta cobertura) o algún grupo externo (FASTA en coordenadas de la
referencia, p. ej. chimpancé proyectado a hg19) lleva un alelo distinto de REF.
Los especímenes de baja cobertura (Sima) NO intervienen en la selección: así sus
errores de secuenciación no crean sitios, solo pueden caer sobre sitios ya elegidos.

    python3 cloud/09_site_panel.py --ref hg19.fa \
        --vcf AltaiNea.vcf.gz Vindija.vcf.gz Denisova.vcf.gz ... \
        --outgroup-fasta panTro4_en_hg19.fa \
        [--chroms chr1 ... chr22] [--mask mapabilidad.bed] [--min-gq 20] \
        [--region chr1:100001-200000] --out pairwise/windows/<id>/panel.tsv.gz

Los VCF tienen que estar comprimidos con bgzip e indexados (tabix).
Para acelerar con VCF de todos los sitios: bcftools view -i 'GT="alt"' antes de pasarlos aquí.
"""
from __future__ import annotations

import argparse
import bisect
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pysam

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
from common import merge_regions, read_bed  # noqa: E402
from pairwise_io import TRANSVERSION_PAIRS, Panel  # noqa: E402


def vcf_alts(vcf: Path, chrom: str, start0: int, end0: int, min_gq: int, stats: dict) -> dict[int, set]:
    out = defaultdict(set)
    with pysam.VariantFile(str(vcf)) as vf:
        if chrom not in vf.header.contigs:
            return out
        for rec in vf.fetch(chrom, start0, end0):
            if len(rec.ref) != 1 or not rec.alts:
                continue
            for s in rec.samples.values():
                gt = s.get("GT")
                if not gt or None in gt:
                    continue
                gq = s.get("GQ") if "GQ" in rec.format else None
                if gq is not None and gq < min_gq:
                    stats["low_gq"] += 1
                    continue
                for i in gt:
                    if i and len(rec.alleles[i]) == 1:
                        out[rec.pos].add(rec.alleles[i].upper())
            if rec.pos in out:
                out[rec.pos].add("REF=" + rec.ref.upper())
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref", required=True, type=Path)
    ap.add_argument("--vcf", nargs="*", default=[], type=Path)
    ap.add_argument("--outgroup-fasta", nargs="*", default=[], type=Path)
    ap.add_argument("--chroms", nargs="*")
    ap.add_argument("--region", help="solo esta ventana, chrom:inicio-fin (1-based, inclusivo); lo usa run_vectors.sh")
    ap.add_argument("--mask", type=Path, help="BED de regiones permitidas (p. ej. mapabilidad alta)")
    ap.add_argument("--min-gq", type=int, default=20)
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args()

    ref = pysam.FastaFile(str(a.ref))
    if a.region:
        chrom, _, span = a.region.rpartition(":")
        b, e = (int(x) for x in span.replace(",", "").split("-"))
        regions = [(chrom, b - 1, e)]
    else:
        regions = [(c, 0, ref.get_reference_length(c)) for c in (a.chroms or ref.references)]
    mask = defaultdict(list)
    if a.mask:
        for r in merge_regions(read_bed(a.mask)):  # fusionar: la búsqueda binaria lo exige
            mask[r.chrom].append((r.start, r.end))
    stats = defaultdict(int)
    sites = []
    for chrom, start0, end0 in regions:
        end0 = min(end0, ref.get_reference_length(chrom))
        refseq = np.frombuffer(ref.fetch(chrom, start0, end0).upper().encode(), dtype=np.uint8)
        stats["ref_N"] += int(np.sum(refseq == ord("N")))
        cand = defaultdict(set)
        for vcf in a.vcf:
            for pos, alleles in vcf_alts(vcf, chrom, start0, end0, a.min_gq, stats).items():
                if not start0 < pos <= end0:
                    continue
                vref = next((x[4:] for x in alleles if x.startswith("REF=")), None)
                if vref and vref != chr(refseq[pos - 1 - start0]):
                    stats["vcf_ref_mismatch"] += 1   # el VCF no está en el ensamblaje de --ref
                    continue
                cand[pos] |= {x for x in alleles if not x.startswith("REF=")}
        for fa in a.outgroup_fasta:
            with pysam.FastaFile(str(fa)) as og:
                if chrom not in og.references:
                    continue
                o = np.frombuffer(og.fetch(chrom, start0, end0).upper().encode(), dtype=np.uint8)
                n = min(len(o), len(refseq))
                valid = np.isin(o[:n], np.frombuffer(b"ACGT", np.uint8)) & (o[:n] != refseq[:n])
                for i in np.nonzero(valid)[0]:
                    cand[start0 + int(i) + 1].add(chr(o[i]))
        ivs = sorted(mask.get(chrom, []))
        starts = [s for s, _ in ivs]

        def in_mask(pos):  # BED 0-based semiabierto: pos (1-based) dentro si s < pos <= e
            i = bisect.bisect_left(starts, pos) - 1
            return i >= 0 and pos <= ivs[i][1]
        for pos in sorted(cand):
            r = chr(refseq[pos - 1 - start0])
            alts = cand[pos] - {r}
            if r not in "ACGT" or len(alts) != 1:
                stats["not_biallelic"] += 1
                continue
            alt = alts.pop()
            if frozenset((r, alt)) not in TRANSVERSION_PAIRS:
                stats["transition"] += 1
                continue
            if a.mask and not in_mask(pos):
                stats["masked"] += 1
                continue
            sites.append((chrom, pos, r, alt))
        print(f"  {chrom}: {sum(1 for s in sites if s[0] == chrom):,} sitios", file=sys.stderr)
    md5 = Panel.write(a.out, sites)
    print(f"Panel: {len(sites):,} transversiones bialélicas -> {a.out} (md5 {md5})")
    print("Descartados (ref_N = bases N de la referencia en la región): " + ", ".join(f"{k}={v:,}" for k, v in sorted(stats.items())))
    if stats["vcf_ref_mismatch"]:
        print("AVISO: hay sitios cuyo REF del VCF no coincide con --ref; revisa el ensamblaje", file=sys.stderr)


if __name__ == "__main__":
    main()
