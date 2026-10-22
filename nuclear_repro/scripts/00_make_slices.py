#!/usr/bin/env python3
"""Paso 0 (se ejecuta UNA vez, donde están los BAM completos): recorta BAM y referencia.

Conserva solo las lecturas que caen en (a) cada posición del catálogo ± WINDOW y
(b) las regiones CDS. Las coordenadas y la cabecera @SQ no cambian, así que los
recortes se analizan exactamente igual que los BAM completos.

    python3 scripts/00_make_slices.py \
        --ref /ruta/hg19.fa \
        --bam /ruta/Sima_FemurXIII.hg19.bam /ruta/Sima_molar.hg19.bam ... \
        [--out slices]

Escribe slices/<muestra>.bam(.bai), slices/regions.fa(.fai) y slices/MANIFEST.tsv
(origen, md5 del BAM completo y del recorte, lecturas). Ese directorio es lo que se versiona.
"""
from __future__ import annotations

import argparse
import datetime
import sys
from pathlib import Path

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (Region, bam_build, die, load_config, md5, merge_regions,  # noqa: E402
                    read_bed, read_positions, sample_name)

ROOT = Path(__file__).resolve().parent.parent


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref", required=True, type=Path, help="FASTA completo de la referencia (con .fai)")
    ap.add_argument("--bam", required=True, nargs="+", type=Path, help="BAM completos (indexados)")
    ap.add_argument("--catalog", type=Path, default=ROOT / "catalog")
    ap.add_argument("--out", type=Path, default=ROOT / "slices")
    ap.add_argument("--config", type=Path, default=ROOT / "config.env")
    a = ap.parse_args()

    cfg = load_config(a.config)
    window = int(cfg["WINDOW"])
    positions = read_positions(a.catalog / "positions.tsv")
    cds = read_bed(a.catalog / "cds_regions.bed")

    regions = [Region(p.chrom, max(0, p.pos - 1 - window), p.pos + window) for p in positions] + cds
    merged = merge_regions(regions)
    span = sum(r.end - r.start for r in merged)
    print(f"{len(positions)} posiciones ± {window} pb + {len(cds)} regiones CDS "
          f"-> {len(merged)} regiones fusionadas, {span:,} pb")

    a.out.mkdir(parents=True, exist_ok=True)
    manifest = [("file", "sample", "build", "source_path", "source_md5", "reads", "md5")]

    # Referencia recortada, cabeceras chrom:inicio-fin (1-based) como samtools faidx
    fa = pysam.FastaFile(str(a.ref))
    with open(a.out / "regions.fa", "w") as out:
        for r in merged:
            if r.chrom not in fa.references:
                die(f"{r.chrom} no está en {a.ref}")
            end = min(r.end, fa.get_reference_length(r.chrom))
            seq = fa.fetch(r.chrom, r.start, end)
            out.write(f">{r.chrom}:{r.start + 1}-{end}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i:i + 60] + "\n")
    pysam.faidx(str(a.out / "regions.fa"))
    manifest.append(("regions.fa", "-", cfg["BUILD"], str(a.ref), "-", "-", md5(a.out / "regions.fa")))

    for bam_path in a.bam:
        name = sample_name(bam_path)
        with pysam.AlignmentFile(str(bam_path)) as bam:
            build = bam_build(bam)
            if build != cfg["BUILD"]:
                die(f"{bam_path} está en {build} y BUILD={cfg['BUILD']}")
            hdr = bam.header.to_dict()
            hdr.setdefault("PG", []).append({
                "ID": "nuclear_repro_slice", "PN": "00_make_slices.py",
                "DS": f"recorte: {len(merged)} regiones, {span} pb; {datetime.date.today()}"})
            tmp = a.out / f"{name}.unsorted.bam"
            n = 0
            seen = set()
            with pysam.AlignmentFile(str(tmp), "wb", header=hdr) as out:
                for r in merged:
                    if r.chrom not in bam.references:
                        continue
                    for read in bam.fetch(r.chrom, r.start, r.end):
                        # una lectura puede tocar dos regiones fusionadas contiguas
                        key = (read.query_name, read.flag, read.reference_start)
                        if key in seen:
                            continue
                        seen.add(key)
                        out.write(read)
                        n += 1
        dest = a.out / f"{name}.bam"
        pysam.sort("-o", str(dest), str(tmp))
        tmp.unlink()
        pysam.index(str(dest))
        manifest.append((dest.name, name, build, str(bam_path), md5(bam_path), str(n), md5(dest)))
        print(f"  {name}: {n:,} lecturas -> {dest} ({dest.stat().st_size / 1e6:.2f} MB)")

    with open(a.out / "MANIFEST.tsv", "w") as fh:
        for row in manifest:
            fh.write("\t".join(row) + "\n")
    print(f"Hecho. Revisa {a.out}/MANIFEST.tsv y versiona el directorio {a.out.name}/.")


if __name__ == "__main__":
    main()
