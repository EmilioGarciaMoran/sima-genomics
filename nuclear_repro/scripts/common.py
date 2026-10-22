"""Utilidades compartidas del pipeline nuclear reproducible.

Coordenadas: el catálogo usa posiciones 1-based (como un VCF); los BED son 0-based
semiabiertos (como siempre). Todo lo demás se convierte en el momento de leerlo.
"""
from __future__ import annotations

import csv
import hashlib
import os
import sys
from dataclasses import dataclass
from pathlib import Path

# Longitud de chr1 en cada ensamblaje: basta para saber contra qué referencia se alineó un BAM.
CHR1_LENGTH = {249250621: "hg19", 248956422: "hg38",
               123456: "fixture"}  # genoma sintético de tests/make_fixture.py
TRANSITIONS = {frozenset("AG"), frozenset("CT")}
BASES = "ACGT"


def die(msg: str) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


def load_config(path: Path) -> dict:
    """Lee config.env (KEY=valor). Las variables de entorno con el mismo nombre mandan."""
    cfg = {}
    for line in Path(path).read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if not line:
            continue
        key, _, val = line.partition("=")
        cfg[key.strip()] = val.strip().strip('"')
    for key in list(cfg):
        if key in os.environ:
            cfg[key] = os.environ[key]
    return cfg


@dataclass(frozen=True)
class Position:
    id: str
    chrom: str
    pos: int  # 1-based
    ref: str
    anc: str
    der: str
    gene: str
    panel: str
    source: str

    @property
    def is_transition(self) -> bool:
        return frozenset((self.anc, self.der)) in TRANSITIONS


POSITION_COLUMNS = ["id", "chrom", "pos", "ref", "anc", "der", "gene", "panel", "source"]


def read_positions(path: Path) -> list[Position]:
    rows = []
    with open(path) as fh:
        reader = csv.DictReader((l for l in fh if not l.startswith("##")), delimiter="\t")
        missing = set(POSITION_COLUMNS) - set(reader.fieldnames or [])
        if missing:
            die(f"{path}: faltan columnas {sorted(missing)}")
        for r in reader:
            for k in ("ref", "anc", "der"):
                r[k] = r[k].upper()
                if r[k] not in BASES:
                    die(f"{path}: {r['id']} tiene {k}={r[k]!r}; solo se admiten SNV (A/C/G/T)")
            rows.append(Position(**{k: (int(r[k]) if k == "pos" else r[k]) for k in POSITION_COLUMNS}))
    if not rows:
        die(f"{path} no tiene posiciones")
    return rows


@dataclass(frozen=True)
class Region:
    chrom: str
    start: int  # 0-based
    end: int    # semiabierto
    name: str = ""


def read_bed(path: Path) -> list[Region]:
    out = []
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        f = line.split("\t")
        out.append(Region(f[0], int(f[1]), int(f[2]), f[3] if len(f) > 3 else ""))
    return out


def merge_regions(regions: list[Region]) -> list[Region]:
    merged: list[Region] = []
    for r in sorted(regions, key=lambda r: (r.chrom, r.start)):
        if merged and merged[-1].chrom == r.chrom and r.start <= merged[-1].end:
            m = merged[-1]
            merged[-1] = Region(m.chrom, m.start, max(m.end, r.end))
        else:
            merged.append(Region(r.chrom, r.start, r.end))
    return merged


class RefSlices:
    """Referencia recortada: FASTA con cabeceras 'chrom:inicio-fin' (1-based, como samtools faidx).

    Permite consultar la base de referencia en coordenadas originales sin
    distribuir el genoma completo.
    """

    def __init__(self, fasta: Path):
        self.slices: dict[str, list[tuple[int, str]]] = {}
        name, seq = None, []
        for line in Path(fasta).read_text().splitlines():
            if line.startswith(">"):
                self._add(name, seq)
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip().upper())
        self._add(name, seq)

    def _add(self, name, seq):
        if name is None:
            return
        chrom, _, span = name.rpartition(":")
        start = int(span.split("-")[0])
        self.slices.setdefault(chrom, []).append((start, "".join(seq)))

    def base(self, chrom: str, pos1: int) -> str | None:
        for start, seq in self.slices.get(chrom, []):
            if start <= pos1 < start + len(seq):
                return seq[pos1 - start]
        return None

    def covers(self, chrom: str, start0: int, end0: int) -> bool:
        return all(self.base(chrom, p + 1) is not None for p in (start0, end0 - 1))


def bam_build(bam) -> str:
    for sq in bam.header.to_dict().get("SQ", []):
        if sq["SN"] in ("chr1", "1"):
            return CHR1_LENGTH.get(sq["LN"], f"desconocido(chr1={sq['LN']})")
    return "desconocido(sin chr1)"


def sample_name(bam_path: Path) -> str:
    return Path(bam_path).name.split(".")[0]


def stable_rng_index(n: int, *key) -> int:
    """Índice reproducible en [0, n): no depende del orden de ejecución ni de PYTHONHASHSEED."""
    h = hashlib.sha256("|".join(map(str, key)).encode()).digest()
    return int.from_bytes(h[:8], "big") % n


def md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def pileup_bases(bam, chrom: str, pos0: int, min_mapq: int, min_bq: int):
    """Bases que cubren una posición: lista de (base, calidad, es_reverse, dist_al_extremo)."""
    out = []
    for col in bam.pileup(chrom, pos0, pos0 + 1, truncate=True, min_mapping_quality=min_mapq,
                          min_base_quality=min_bq, ignore_orphans=False, stepper="nofilter"):
        for pr in col.pileups:
            if pr.is_del or pr.is_refskip or pr.query_position is None:
                continue
            aln = pr.alignment
            if aln.is_unmapped or aln.is_secondary or aln.is_duplicate or aln.is_qcfail:
                continue
            qp = pr.query_position
            b = aln.query_sequence[qp].upper()
            if b not in BASES:
                continue
            dist = min(qp, aln.query_length - 1 - qp)
            out.append((b, aln.query_qualities[qp], aln.is_reverse, dist))
    return out
