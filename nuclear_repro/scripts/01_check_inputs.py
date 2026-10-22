#!/usr/bin/env python3
"""Paso 1: comprueba que los datos de entrada son los que dicen ser antes de analizar nada.

Aborta si:
  - algún recorte no coincide con el md5 de slices/MANIFEST.tsv;
  - algún BAM está en otro ensamblaje que BUILD (longitud de chr1 en la cabecera);
  - el alelo REF de alguna posición del catálogo no coincide con la referencia recortada
    (síntoma típico de mezclar coordenadas hg19 y hg38);
  - alguna posición o región CDS cae fuera de lo recortado.
"""
import sys
from pathlib import Path

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import RefSlices, bam_build, die, load_config, md5, read_bed, read_positions  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent


def main(slices: Path, catalog: Path, cfg_path: Path) -> None:
    cfg = load_config(cfg_path)
    problems = []

    manifest = slices / "MANIFEST.tsv"
    if not manifest.exists():
        die(f"no existe {manifest} (se genera con 00_make_slices.py)")
    rows = [l.split("\t") for l in manifest.read_text().splitlines()[1:] if l.strip()]
    bams = []
    for f, sample, build, _src, _src_md5, _reads, digest in rows:
        p = slices / f
        if not p.exists():
            problems.append(f"falta {p}")
            continue
        if md5(p) != digest.strip():
            problems.append(f"{f}: md5 distinto del de MANIFEST.tsv (¿se ha modificado el recorte?)")
        if f.endswith(".bam"):
            bams.append(p)
    print(f"MANIFEST: {len(rows)} ficheros, {len(bams)} BAM")

    for p in bams:
        with pysam.AlignmentFile(str(p)) as bam:
            b = bam_build(bam)
        if b != cfg["BUILD"]:
            problems.append(f"{p.name}: ensamblaje {b}, BUILD={cfg['BUILD']}")
    print(f"Ensamblaje esperado: {cfg['BUILD']}")

    ref = RefSlices(slices / "regions.fa")
    positions = read_positions(catalog / "positions.tsv")
    mism = []
    for pos in positions:
        base = ref.base(pos.chrom, pos.pos)
        if base is None:
            problems.append(f"{pos.id} ({pos.chrom}:{pos.pos}) no está en regions.fa")
        elif base != pos.ref:
            mism.append(f"{pos.id} {pos.chrom}:{pos.pos} catálogo={pos.ref} referencia={base}")
        if pos.ref not in (pos.anc, pos.der):
            print(f"  aviso: {pos.id} REF={pos.ref} no es ni el ancestral ({pos.anc}) ni el derivado ({pos.der})")
    if mism:
        problems.append(f"{len(mism)}/{len(positions)} posiciones con REF distinto de la referencia "
                        f"(¿catálogo y BAM en ensamblajes distintos?):\n    " + "\n    ".join(mism[:10]))
    print(f"Catálogo: {len(positions)} posiciones, REF verificado en {len(positions) - len(mism)}")

    cds = read_bed(catalog / "cds_regions.bed")
    outside = [f"{r.name} {r.chrom}:{r.start}-{r.end}" for r in cds if not ref.covers(r.chrom, r.start, r.end)]
    if outside:
        problems.append(f"{len(outside)} regiones CDS fuera de regions.fa: " + ", ".join(outside[:5]))
    print(f"CDS: {len(cds)} regiones, {sum(r.end - r.start for r in cds):,} pb")

    if problems:
        die("entradas no válidas:\n  - " + "\n  - ".join(problems))
    print("OK: entradas verificadas")


if __name__ == "__main__":
    slices = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "slices"
    catalog = Path(sys.argv[2]) if len(sys.argv) > 2 else ROOT / "catalog"
    main(slices, catalog, ROOT / "config.env")
