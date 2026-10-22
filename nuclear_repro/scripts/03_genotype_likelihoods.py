#!/usr/bin/env python3
"""Paso 3: verosimilitudes genotípicas sobre las regiones CDS, sitio a sitio.

Modelo (el de GATK, equivalente a ANGSD -GL 2): para el genotipo a1a2 y una base b
con calidad Q (error e = 10^(-Q/10)),
    P(b | a1a2) = ½ P(b|a1) + ½ P(b|a2),   P(b|a) = 1-e si b = a, e/3 si no.
Mayor/menor: el par de alelos que maximiza la verosimilitud (aproximación a
ANGSD -doMajorMinor 1). Se dan log10 GL normalizadas (máximo = 0) de MM, Mm y mm.

Está en Python puro para que el flujo se pueda leer y probar sin ANGSD; el paso 3b
repite el análisis con ANGSD y los parámetros exactos del manuscrito si está instalado.

Salidas: results/gl_sites.tsv (sitios con ≥1 lectura) y results/gl_summary.tsv (por gen).
"""
import itertools
import math
import sys
from collections import defaultdict
from pathlib import Path

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import BASES, TRANSITIONS, RefSlices, load_config, read_bed, sample_name  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
PAIRS = list(itertools.combinations(BASES, 2))


def log10_gl(obs, a1, a2):
    total = 0.0
    for b, q in obs:
        e = 10 ** (-max(q, 1) / 10)
        p1 = 1 - e if b == a1 else e / 3
        p2 = 1 - e if b == a2 else e / 3
        total += math.log10(0.5 * p1 + 0.5 * p2)
    return total


def site_gl(obs):
    """Devuelve (mayor, menor, [GL_MM, GL_Mm, GL_mm] normalizadas)."""
    best = None
    for x, y in PAIRS:
        gls = [log10_gl(obs, x, x), log10_gl(obs, x, y), log10_gl(obs, y, y)]
        m = max(gls)
        if best is None or m > best[0]:
            best = (m, x, y, gls)
    _, x, y, gls = best
    # mayor = el alelo con más apoyo ponderado por calidad (desempate alfabético, determinista)
    w = defaultdict(float)
    for b, q in obs:
        w[b] += q
    if (w[y], y) > (w[x], x):
        x, y, gls = y, x, gls[::-1]
    m = max(gls)
    return x, y, [g - m for g in gls]


def collect(bams, region, min_mapq, min_bq, trim):
    """{pos0: [(base, qual, dist_extremo)]} para la región, juntando todos los BAM del grupo."""
    sites = defaultdict(list)
    for bam in bams:
        if region.chrom not in bam.references:
            continue
        for col in bam.pileup(region.chrom, region.start, region.end, truncate=True,
                              min_mapping_quality=min_mapq, min_base_quality=min_bq,
                              ignore_orphans=False, stepper="nofilter"):
            for pr in col.pileups:
                if pr.is_del or pr.is_refskip or pr.query_position is None:
                    continue
                aln = pr.alignment
                if aln.is_secondary or aln.is_duplicate or aln.is_qcfail:
                    continue
                qp = pr.query_position
                b = aln.query_sequence[qp].upper()
                if b in BASES:
                    sites[col.reference_pos].append((b, aln.query_qualities[qp],
                                                     min(qp, aln.query_length - 1 - qp)))
    return sites


def main(slices: Path, catalog: Path, out: Path) -> None:
    cfg = load_config(ROOT / "config.env")
    mq, bq, trim = int(cfg["MIN_MAPQ_GL"]), int(cfg["MIN_BQ"]), int(cfg["END_TRIM"])
    regions = read_bed(catalog / "cds_regions.bed")
    ref = RefSlices(slices / "regions.fa")
    bams = {sample_name(p): pysam.AlignmentFile(str(p)) for p in sorted(slices.glob("*.bam"))}
    out.mkdir(parents=True, exist_ok=True)

    cols = ["sample", "gene", "chrom", "pos", "ref", "depth", "A", "C", "G", "T", "major", "minor",
            "transition", "GL_MM", "GL_Mm", "GL_mm", "best_gt", "gt_bases", "nonref_reads", "nonref_reads_trimmed"]
    summ = []
    with open(out / "gl_sites.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for sample in list(bams) + ["POOL"]:
            group = list(bams.values()) if sample == "POOL" else [bams[sample]]
            per_gene = defaultdict(lambda: dict(bp=0, cov=0, depth=0, nonref=0, nonref_tv=0, nonref_tv_trim=0))
            for reg in regions:
                g = per_gene[reg.name]
                g["bp"] += reg.end - reg.start
                sites = collect(group, reg, mq, bq, trim)
                for pos0 in sorted(sites):
                    obs = sites[pos0]
                    r = ref.base(reg.chrom, pos0 + 1) or "N"
                    major, minor, gls = site_gl([(b, q) for b, q, _ in obs])
                    counts = {x: sum(1 for b, _, _ in obs if b == x) for x in BASES}
                    nonref = [o for o in obs if o[0] != r]
                    nonref_trim = [o for o in nonref if o[2] >= trim]
                    best = ["MM", "Mm", "mm"][gls.index(0.0)]
                    gt_bases = {"MM": major * 2, "Mm": "".join(sorted(major + minor)), "mm": minor * 2}[best]
                    # transición/transversión respecto a la referencia, según la base no-REF observada
                    alt = max((x for x in BASES if x != r), key=lambda x: (counts[x], x))
                    is_ts = frozenset((r, alt)) in TRANSITIONS
                    fh.write("\t".join(map(str, [
                        sample, reg.name, reg.chrom, pos0 + 1, r, len(obs),
                        counts["A"], counts["C"], counts["G"], counts["T"], major, minor,
                        int(frozenset((major, minor)) in TRANSITIONS),
                        f"{gls[0]:.3f}", f"{gls[1]:.3f}", f"{gls[2]:.3f}", best, gt_bases,
                        len(nonref), len(nonref_trim)])) + "\n")
                    g["cov"] += 1
                    g["depth"] += len(obs)
                    if nonref:
                        g["nonref"] += 1
                        if not is_ts:
                            g["nonref_tv"] += 1
                            if nonref_trim:
                                g["nonref_tv_trim"] += 1
            for gene, g in sorted(per_gene.items()):
                summ.append((sample, gene, g))
            tot = {k: sum(g[k] for g in per_gene.values()) for k in ("bp", "cov", "nonref", "nonref_tv")}
            print(f"  {sample}: {tot['cov']:,}/{tot['bp']:,} pb CDS cubiertos; sitios con base no-REF: "
                  f"{tot['nonref']} (transversiones: {tot['nonref_tv']})")

    with open(out / "gl_summary.tsv", "w") as fh:
        fh.write("sample\tgene\tcds_bp\tbp_covered\tfrac_covered\tmean_depth_covered\t"
                 "sites_nonref\tsites_nonref_transversion\tsites_nonref_transversion_trimmed\n")
        for sample, gene, g in summ:
            frac = g["cov"] / g["bp"] if g["bp"] else 0
            md = g["depth"] / g["cov"] if g["cov"] else 0
            fh.write(f"{sample}\t{gene}\t{g['bp']}\t{g['cov']}\t{frac:.4f}\t{md:.2f}\t"
                     f"{g['nonref']}\t{g['nonref_tv']}\t{g['nonref_tv_trim']}\n")


if __name__ == "__main__":
    a = sys.argv[1:] + [None] * 3
    main(Path(a[0] or ROOT / "slices"), Path(a[1] or ROOT / "catalog"), Path(a[2] or ROOT / "results"))
