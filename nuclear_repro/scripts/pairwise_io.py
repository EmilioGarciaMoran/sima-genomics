"""Formato del panel de sitios informativos y de los vectores por espécimen.

Panel (pairwise/panel.tsv.gz): una fila por sitio, ordenadas por cromosoma y posición
    chrom  pos(1-based)  ref  alt          solo SNV bialélicos de transversión

Vector (pairwise/vectors/<espécimen>.vec.gz): el alelo del espécimen en cada sitio del
panel, en el mismo orden, un carácter por sitio y una línea por cromosoma:
    #specimen=...  #panel_md5=...  #mode=...  (cabeceras de procedencia)
    chr1<TAB>AC.NNA...
con '0' = REF, '1' = ALT y 'N' = sin dato (o un alelo que no es ni REF ni ALT).

El md5 del panel se calcula sobre el texto descomprimido: los gz se escriben con
mtime=0 para que el mismo contenido dé el mismo fichero byte a byte.
"""
from __future__ import annotations

import gzip
import hashlib
import io
from collections import OrderedDict
from pathlib import Path

import numpy as np

CODE = {"0": 1, "1": 2, "N": 0}  # en memoria: 0 = sin dato, 1 = REF, 2 = ALT
TRANSVERSION_PAIRS = {frozenset(p) for p in ("AC", "AT", "GC", "GT")}


def write_gz_text(path: Path, text: str) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0, filename="") as gz:
        gz.write(text.encode())


def read_gz_text(path: Path) -> str:
    with gzip.open(path, "rt") as fh:
        return fh.read()


class Panel:
    def __init__(self, chroms: "OrderedDict[str, np.ndarray]", ref: dict, alt: dict, md5: str):
        self.chroms, self.ref, self.alt, self.md5 = chroms, ref, alt, md5

    @property
    def n(self) -> int:
        return sum(len(p) for p in self.chroms.values())

    @classmethod
    def read(cls, path: Path) -> "Panel":
        text = read_gz_text(path)
        md5 = hashlib.md5(text.encode()).hexdigest()
        pos, ref, alt = OrderedDict(), {}, {}
        for line in text.splitlines():
            if not line or line.startswith(("#", "chrom\t")):
                continue
            c, p, r, a = line.split("\t")[:4]
            pos.setdefault(c, []).append(int(p))
            ref.setdefault(c, []).append(r)
            alt.setdefault(c, []).append(a)
        chroms = OrderedDict((c, np.array(v, dtype=np.int64)) for c, v in pos.items())
        for c, arr in chroms.items():
            if np.any(np.diff(arr) <= 0):
                raise ValueError(f"{path}: posiciones de {c} no ordenadas o repetidas")
        return cls(chroms, ref, alt, md5)

    @staticmethod
    def write(path: Path, sites) -> str:
        """sites: iterable de (chrom, pos, ref, alt) ya ordenado. Devuelve el md5."""
        buf = io.StringIO()
        buf.write("chrom\tpos\tref\talt\n")
        for c, p, r, a in sites:
            buf.write(f"{c}\t{p}\t{r}\t{a}\n")
        text = buf.getvalue()
        write_gz_text(path, text)
        return hashlib.md5(text.encode()).hexdigest()


def write_vector(path: Path, panel: Panel, alleles: dict, meta: dict) -> None:
    """alleles: {chrom: lista de '0'/'1'/'N' en el orden del panel}."""
    lines = [f"#{k}={v}" for k, v in {"panel_md5": panel.md5, **meta}.items()]
    for c in panel.chroms:
        lines.append(f"{c}\t{''.join(alleles[c])}")
    write_gz_text(path, "\n".join(lines) + "\n")


def read_vector(path: Path, panel: Panel) -> tuple[dict, np.ndarray]:
    meta, parts = {}, {}
    for line in read_gz_text(path).splitlines():
        if line.startswith("#"):
            k, _, v = line[1:].partition("=")
            meta[k] = v
        elif line:
            c, s = line.split("\t")
            parts[c] = s
    if meta.get("panel_md5") != panel.md5:
        raise ValueError(f"{path}: hecho con otro panel (panel_md5 {meta.get('panel_md5')} ≠ {panel.md5})")
    out = []
    for c, arr in panel.chroms.items():
        s = parts.get(c, "")
        if len(s) != len(arr):
            raise ValueError(f"{path}: {c} tiene {len(s)} sitios y el panel {len(arr)}")
        out.append(np.frombuffer(s.translate(str.maketrans("01N", "\x01\x02\x00")).encode("latin1"),
                                 dtype=np.uint8))
    return meta, np.concatenate(out) if out else np.zeros(0, np.uint8)


def encode_allele(base: str | None, ref: str, alt: str) -> str:
    return "0" if base == ref else "1" if base == alt else "N"
