#!/usr/bin/env python3
"""Paso 2: genotipado directo de las posiciones del catálogo en cada BAM y en el conjunto.

Para cada posición y muestra (más "POOL", todas las lecturas juntas):
  - recuento de A/C/G/T tras filtros (MIN_MAPQ_DIRECT, MIN_BQ);
  - llamada pseudohaploide: una lectura elegida al azar con semilla fija, como es
    habitual a baja cobertura; se clasifica como ANC, DER, OTRO o NA (sin lectura);
  - la misma llamada descartando las END_TRIM bases de cada extremo (daño post mortem);
  - lecturas en la ventana ± WINDOW (para cuantificar "hay lecturas cerca pero no encima").
Las posiciones de transición (C/T, G/A) se marcan: una llamada DER ahí puede ser daño.

Salidas: results/direct_genotypes.tsv y results/direct_summary.tsv
"""
import sys
from pathlib import Path

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import load_config, pileup_bases, read_positions, sample_name, stable_rng_index  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent


def classify(base, pos):
    if base is None:
        return "NA"
    return "ANC" if base == pos.anc else "DER" if base == pos.der else "OTRO"


def call(bases, pos, sample, seed, tag):
    if not bases:
        return None
    i = stable_rng_index(len(bases), seed, sample, pos.chrom, pos.pos, tag)
    return sorted(bases)[i][0]  # orden fijo antes de elegir: reproducible entre versiones de htslib


def reads_in_window(bams, chrom, pos0, window, min_mapq):
    n = 0
    for bam in bams:
        if chrom not in bam.references:
            continue
        for r in bam.fetch(chrom, max(0, pos0 - window), pos0 + window + 1):
            if r.mapping_quality >= min_mapq and not (r.is_unmapped or r.is_secondary or r.is_duplicate):
                n += 1
    return n


def main(slices: Path, catalog: Path, out: Path) -> None:
    cfg = load_config(ROOT / "config.env")
    mq, bq, win, trim, seed = (int(cfg[k]) for k in ("MIN_MAPQ_DIRECT", "MIN_BQ", "WINDOW", "END_TRIM", "SEED"))
    positions = read_positions(catalog / "positions.tsv")
    bam_paths = sorted(slices.glob("*.bam"))
    bams = {sample_name(p): pysam.AlignmentFile(str(p)) for p in bam_paths}

    out.mkdir(parents=True, exist_ok=True)
    cols = ["sample", "id", "gene", "panel", "chrom", "pos", "anc", "der", "transition",
            "depth", "A", "C", "G", "T", "call", "state", "depth_trimmed", "call_trimmed",
            "state_trimmed", "reads_window"]
    summary = {}
    with open(out / "direct_genotypes.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for sample in list(bams) + ["POOL"]:
            group = list(bams.values()) if sample == "POOL" else [bams[sample]]
            s = summary.setdefault(sample, dict(n=0, cov1=0, cov2=0, near=0, anc=0, der=0, other=0, der_tv=0))
            for pos in positions:
                bases = []
                for bam in group:
                    if pos.chrom in bam.references:
                        bases += pileup_bases(bam, pos.chrom, pos.pos - 1, mq, bq)
                trimmed = [b for b in bases if b[3] >= trim]
                counts = {x: sum(1 for b in bases if b[0] == x) for x in "ACGT"}
                c = call(bases, pos, sample, seed, "all")
                ct = call(trimmed, pos, sample, seed, "trim")
                near = reads_in_window(group, pos.chrom, pos.pos - 1, win, mq)
                st = classify(c, pos)
                fh.write("\t".join(map(str, [
                    sample, pos.id, pos.gene, pos.panel, pos.chrom, pos.pos, pos.anc, pos.der,
                    int(pos.is_transition), len(bases), counts["A"], counts["C"], counts["G"], counts["T"],
                    c or ".", st, len(trimmed), ct or ".", classify(ct, pos), near])) + "\n")
                s["n"] += 1
                s["cov1"] += len(bases) >= 1
                s["cov2"] += len(bases) >= 2
                s["near"] += near > 0
                s["anc"] += st == "ANC"
                s["der"] += st == "DER"
                s["other"] += st == "OTRO"
                s["der_tv"] += st == "DER" and not pos.is_transition

    with open(out / "direct_summary.tsv", "w") as fh:
        fh.write("sample\tpositions\tcovered_ge1\tcovered_ge2\treads_within_window\tANC\tDER\tOTRO\tDER_transversion\n")
        for sample, s in summary.items():
            fh.write("\t".join(map(str, [sample, s["n"], s["cov1"], s["cov2"], s["near"],
                                         s["anc"], s["der"], s["other"], s["der_tv"]])) + "\n")
            print(f"  {sample}: {s['cov1']}/{s['n']} posiciones con ≥1 lectura, "
                  f"{s['near']}/{s['n']} con lecturas a ±{win} pb; ANC={s['anc']} DER={s['der']}")


if __name__ == "__main__":
    a = sys.argv[1:] + [None] * 3
    main(Path(a[0] or ROOT / "slices"), Path(a[1] or ROOT / "catalog"), Path(a[2] or ROOT / "results"))
