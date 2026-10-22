#!/usr/bin/env python3
"""NUBE — paso 10: el vector de alelos de un espécimen sobre el panel de sitios.

Un solo alelo por sitio para todos los especímenes, también los de alta cobertura,
para que las distancias sean comparables entre sí:
  --bam    baja cobertura (Sima, humanos antiguos): una lectura al azar por sitio
           (pseudohaploide), filtros de MAPQ/calidad de base, opcionalmente sin extremos.
  --vcf    alta cobertura: uno de los dos alelos del genotipo al azar (o el único si es homocigoto).
  --fasta  grupo externo proyectado a la referencia (chimpancé, gorila, orangután).
Un alelo que no es ni REF ni ALT del panel (error, daño, tercer alelo) se guarda como N.
El azar es reproducible: depende de SEED, del espécimen y del sitio.

    python3 cloud/10_specimen_vector.py --panel pairwise/windows/<id>/panel.tsv.gz --name Sima_FemurXIII \
        --bam Sima_FemurXIII.hg19.bam --out pairwise/windows/<id>/vectors
    python3 cloud/10_specimen_vector.py --panel ... --name Vindija33.19 --vcf Vindija.vcf.gz [--sample X]
    python3 cloud/10_specimen_vector.py --panel ... --name Chimp --fasta panTro4_en_hg19.fa
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pysam

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
from common import load_config, stable_rng_index  # noqa: E402
from pairwise_io import Panel, encode_allele, write_vector  # noqa: E402


def from_bam(path, panel, name, cfg):
    mq, bq, trim, seed = (int(cfg[k]) for k in ("PAIR_MIN_MAPQ", "MIN_BQ", "PAIR_END_TRIM", "SEED"))
    out, stats = {}, dict(reads=0, sites_covered=0)
    with pysam.AlignmentFile(str(path)) as bam:
        for chrom, pos in panel.chroms.items():
            bases = [[] for _ in range(len(pos))]
            if chrom in bam.references and len(pos):
                # solo el tramo del panel (una ventana): no se recorre el cromosoma entero
                for r in bam.fetch(chrom, int(pos[0]) - 1, int(pos[-1])):
                    if (r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate
                            or r.is_qcfail or r.mapping_quality < mq):
                        continue
                    # ¿cae algún sitio dentro de la lectura? (posiciones 1-based en el panel)
                    lo = np.searchsorted(pos, r.reference_start + 1)
                    hi = np.searchsorted(pos, r.reference_end, side="right")
                    if lo >= hi:
                        continue
                    stats["reads"] += 1
                    want = {int(pos[i]) - 1: i for i in range(lo, hi)}
                    seq, qual, L = r.query_sequence, r.query_qualities, r.query_length
                    for qp, rp in r.get_aligned_pairs(matches_only=True):
                        i = want.get(rp)
                        if i is None or qual[qp] < bq or min(qp, L - 1 - qp) < trim:
                            continue
                        bases[i].append((seq[qp].upper(), r.query_name))
            alle = []
            for i, b in enumerate(bases):
                if not b:
                    alle.append("N")
                    continue
                stats["sites_covered"] += 1
                b.sort()  # orden fijo antes de elegir
                pick = b[stable_rng_index(len(b), seed, name, chrom, int(pos[i]))][0]
                alle.append(encode_allele(pick, panel.ref[chrom][i], panel.alt[chrom][i]))
            out[chrom] = alle
    return out, stats


def from_vcf(path, panel, name, cfg, sample, min_gq):
    seed = int(cfg["SEED"])
    out, stats = {}, dict(records=0, low_gq=0, het=0)
    with pysam.VariantFile(str(path)) as vf:
        sample = list(vf.header.samples)[0] if sample in (None, "", "-") else sample
        stats["sample"] = sample
        for chrom, pos in panel.chroms.items():
            alle = ["N"] * len(pos)
            if chrom in vf.header.contigs:
                index = {int(p): i for i, p in enumerate(pos)}
                # recorre solo los sitios del panel, por tramos, aprovechando el índice tabix
                for start in range(0, len(pos), 5000):
                    block = pos[start:start + 5000]
                    for rec in vf.fetch(chrom, int(block[0]) - 1, int(block[-1])):
                        i = index.get(rec.pos)
                        if i is None or len(rec.ref) != 1:
                            continue
                        s = rec.samples[sample]
                        gt = s.get("GT")
                        if not gt or None in gt:
                            continue
                        if "GQ" in rec.format and s.get("GQ") is not None and s.get("GQ") < min_gq:
                            stats["low_gq"] += 1
                            continue
                        stats["records"] += 1
                        if len(set(gt)) > 1:
                            stats["het"] += 1
                        k = gt[stable_rng_index(len(gt), seed, name, chrom, rec.pos)]
                        alle[i] = encode_allele(rec.alleles[k].upper(), panel.ref[chrom][i], panel.alt[chrom][i])
            out[chrom] = alle
    return out, stats


def from_fasta(path, panel):
    out, stats = {}, dict(sites_called=0)
    with pysam.FastaFile(str(path)) as fa:
        for chrom, pos in panel.chroms.items():
            if chrom not in fa.references:
                out[chrom] = ["N"] * len(pos)
                continue
            lo = int(pos[0]) - 1 if len(pos) else 0
            seq = fa.fetch(chrom, lo, int(pos[-1])).upper() if len(pos) else ""
            alle = [encode_allele(seq[p - 1 - lo] if p - 1 - lo < len(seq) else None,
                                  panel.ref[chrom][i], panel.alt[chrom][i]) for i, p in enumerate(pos)]
            stats["sites_called"] += sum(x != "N" for x in alle)
            out[chrom] = alle
    return out, stats


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panel", required=True, type=Path)
    ap.add_argument("--name", required=True)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--bam", type=Path)
    src.add_argument("--vcf", type=Path)
    src.add_argument("--fasta", type=Path)
    ap.add_argument("--sample", help="muestra del VCF (por defecto, la primera)")
    ap.add_argument("--min-gq", type=int, default=20)
    ap.add_argument("--out", type=Path, default=ROOT / "pairwise" / "vectors")
    ap.add_argument("--config", type=Path, default=ROOT / "config.env")
    a = ap.parse_args()

    cfg = load_config(a.config)
    panel = Panel.read(a.panel)
    if a.bam:
        alle, stats = from_bam(a.bam, panel, a.name, cfg)
        meta = dict(mode="bam", source=a.bam.name, min_mapq=cfg["PAIR_MIN_MAPQ"], min_bq=cfg["MIN_BQ"],
                    end_trim=cfg["PAIR_END_TRIM"], seed=cfg["SEED"])
    elif a.vcf:
        alle, stats = from_vcf(a.vcf, panel, a.name, cfg, a.sample, a.min_gq)
        meta = dict(mode="vcf", source=a.vcf.name, sample=stats.pop("sample"), min_gq=a.min_gq, seed=cfg["SEED"])
    else:
        alle, stats = from_fasta(a.fasta, panel)
        meta = dict(mode="fasta", source=a.fasta.name)
    called = sum(x != "N" for v in alle.values() for x in v)
    meta = {"specimen": a.name, **meta, "sites": panel.n, "called": called}
    dest = a.out / f"{a.name}.vec.gz"
    write_vector(dest, panel, alle, meta)
    print(f"{a.name}: {called:,}/{panel.n:,} sitios con alelo ({100 * called / max(panel.n, 1):.1f}%) -> {dest}  "
          + " ".join(f"{k}={v:,}" for k, v in stats.items()))


if __name__ == "__main__":
    main()
