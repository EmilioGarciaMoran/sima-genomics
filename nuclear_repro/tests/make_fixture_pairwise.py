#!/usr/bin/env python3
"""Caso sintético para las distancias por pares: datos simulados sobre un árbol conocido.

                ┌── Sima      (BAM ~1x, 2% de error por base, daño C->T en los extremos)
             ┌──┤
          ┌──┤  └── Nea       (VCF diploide, 5% de heterocigotos)
          │  └──── Den        (VCF)
       ┌──┤
       │  └─────── Human      (= la referencia; VCF)
    ┌──┤
    │  └────────── Chimp      (FASTA en coordenadas de la referencia; grupo externo)
    └───────────── Gorilla    (FASTA)

Además 'Sima_true': los alelos verdaderos de Sima sin error, en FASTA, para comprobar
que los estadísticos de Sima (baja cobertura, con error) se parecen a los verdaderos.

Cada sitio lleva una mutación de transversión en una rama elegida con probabilidad
proporcional a su longitud. Lo esperado: Sima comparte alelos derivados con Nea y no
con Den (rama Sima+Nea), así que D(Sima; Nea, Den) >> 0.

    python3 tests/make_fixture_pairwise.py OUTDIR
      -> OUTDIR/data/{ref.fa, Nea.vcf.gz, Den.vcf.gz, Human.vcf.gz, Sima.bam, Chimp.fa, Gorilla.fa, Sima_true.fa}
         OUTDIR/pairwise/specimens.tsv
"""
import random
import sys
from pathlib import Path

import pysam

LEN = {"chr1": 123456, "chr2": 120000}
SPACING = 55
RL = 40
TAXA = ["Sima", "Nea", "Den", "Human", "Chimp", "Gorilla"]
BRANCHES = [  # (descendientes, longitud)
    ({"Sima"}, 1.0), ({"Nea"}, 2.0), ({"Sima", "Nea"}, 1.5), ({"Den"}, 3.0),
    ({"Sima", "Nea", "Den"}, 2.0), ({"Human"}, 3.0), ({"Sima", "Nea", "Den", "Human"}, 6.0),
    ({"Chimp"}, 6.0), ({"Sima", "Nea", "Den", "Human", "Chimp"}, 2.0), ({"Gorilla"}, 8.0)]
PARTNERS = {"A": "CT", "G": "CT", "C": "AG", "T": "AG"}  # transversiones
rng = random.Random(12345)


def write_fa(path, seqs):
    with open(path, "w") as fh:
        for c, s in seqs.items():
            fh.write(f">{c}\n" + "\n".join(s[i:i + 60] for i in range(0, len(s), 60)) + "\n")
    pysam.faidx(str(path))


def write_vcf(path, name, sites, ref, het=0.0):
    hdr = pysam.VariantHeader()
    for c, n in LEN.items():
        hdr.contigs.add(c, length=n)
    hdr.add_meta("FORMAT", items=[("ID", "GT"), ("Number", 1), ("Type", "String"), ("Description", "Genotype")])
    hdr.add_meta("FORMAT", items=[("ID", "GQ"), ("Number", 1), ("Type", "Integer"), ("Description", "GQ")])
    hdr.add_sample(name)
    raw = path.with_suffix("")  # .vcf
    with pysam.VariantFile(str(raw), "w", header=hdr) as vf:
        for (c, p), alle in sites:
            r = ref[c][p - 1]
            a = alle[name]
            alt = a if a != r else rng.choice(PARTNERS[r])
            if rng.random() < het:
                gt = (0, 1)
            else:
                gt = (0, 0) if a == r else (1, 1)
            rec = vf.new_record(contig=c, start=p - 1, stop=p, alleles=(r, alt))
            rec.samples[name]["GT"] = gt
            rec.samples[name]["GQ"] = 50
            vf.write(rec)
    pysam.tabix_index(str(raw), preset="vcf", force=True)


def main(out: Path):
    data = out / "data"
    data.mkdir(parents=True, exist_ok=True)
    base = {c: [rng.choice("ACGT") for _ in range(n)] for c, n in LEN.items()}
    weights = [w for _, w in BRANCHES]
    sites = []
    for c, n in LEN.items():
        for p in range(200, n - 200, SPACING):
            anc = base[c][p - 1]
            der = rng.choice(PARTNERS[anc])
            desc, _ = rng.choices(BRANCHES, weights)[0]
            sites.append(((c, p), {t: der if t in desc else anc for t in TAXA}))
    ref = {c: list(s) for c, s in base.items()}
    for (c, p), alle in sites:            # la referencia es el humano
        ref[c][p - 1] = alle["Human"]
    ref = {c: "".join(s) for c, s in ref.items()}
    write_fa(data / "ref.fa", ref)
    for t in ("Chimp", "Gorilla", "Sima"):
        seq = {c: list(s) for c, s in ref.items()}
        for (c, p), alle in sites:
            seq[c][p - 1] = alle[t]
        write_fa(data / (f"{t}.fa" if t != "Sima" else "Sima_true.fa"), {c: "".join(s) for c, s in seq.items()})
    write_vcf(data / "Nea.vcf.gz", "Nea", sites, ref, het=0.05)
    write_vcf(data / "Den.vcf.gz", "Den", sites, ref)
    write_vcf(data / "Human.vcf.gz", "Human", sites, ref)

    # Sima: lecturas de 40 pb, número de lecturas por sitio ~ Poisson(1)
    header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6"}, "SQ": [{"SN": c, "LN": n} for c, n in LEN.items()]})
    tmp = data / "Sima.unsorted.bam"
    k = 0
    with pysam.AlignmentFile(str(tmp), "wb", header=header) as bam:
        for (c, p), alle in sites:
            nreads = sum(1 for _ in range(3) if rng.random() < 1 / 3)  # Binomial(3, 1/3), media 1
            for _ in range(nreads):
                start = p - 1 - rng.randrange(0, RL)
                seq = list(ref[c][start:start + RL])
                seq[p - 1 - start] = alle["Sima"]
                for i in range(RL):
                    if rng.random() < 0.02:                     # error de secuenciación
                        seq[i] = rng.choice([x for x in "ACGT" if x != seq[i]])
                    if i < 2 and seq[i] == "C" and rng.random() < 0.3:  # daño en el extremo 5'
                        seq[i] = "T"
                a = pysam.AlignedSegment(header)
                a.query_name, a.reference_name, a.reference_start = f"s{k}", c, start
                a.flag, a.mapping_quality, a.cigarstring = 0, 37, f"{RL}M"
                a.query_sequence = "".join(seq)
                a.query_qualities = pysam.qualitystring_to_array("?" * RL)
                bam.write(a)
                k += 1
    pysam.sort("-o", str(data / "Sima.bam"), str(tmp))
    tmp.unlink()
    pysam.index(str(data / "Sima.bam"))

    pw = out / "pairwise"
    (pw / "vectors").mkdir(parents=True, exist_ok=True)
    (pw / "specimens.tsv").write_text(
        "## caso sintético\n"
        "name\tgroup\trole\tmode\tpath\tsample\n"
        "Sima\tSima\ttarget\tbam\tSima.bam\t-\n"
        "Sima_true\tSima\ttarget\tfasta\tSima_true.fa\t-\n"
        "Nea\tNeandertal\tcomparator\tvcf\tNea.vcf.gz\t-\n"
        "Den\tDenisova\tcomparator\tvcf\tDen.vcf.gz\t-\n"
        "Human\tHumano\tcomparator\tvcf\tHuman.vcf.gz\t-\n"
        "Chimp\tSimio\toutgroup\tfasta\tChimp.fa\t-\n"
        "Gorilla\tSimio\tcomparator\tfasta\tGorilla.fa\t-\n")
    print(f"{len(sites)} sitios simulados, {k} lecturas de Sima -> {out}")


if __name__ == "__main__":
    main(Path(sys.argv[1]))
