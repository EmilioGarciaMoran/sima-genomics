#!/usr/bin/env python3
"""Lector de FASTA con acceso aleatorio mediante el índice .fai de samtools.

Evita depender de pysam, que no está disponible en todos los entornos. El .fai lo genera
`samtools faidx` y tiene cinco columnas: nombre, longitud, desplazamiento del primer byte
de secuencia, bases por línea y bytes por línea (bases más el salto de línea).

    ref = Fasta("human_g1k_v37.fasta")
    ref.fetch("21", 9411194, 9411204)   # 0-based, semiabierto, como en Python
"""
from __future__ import annotations

from pathlib import Path


class Fasta:
    def __init__(self, path):
        self.path = Path(path)
        fai = Path(str(path) + ".fai")
        if not fai.exists():
            raise FileNotFoundError(
                f"falta el índice {fai}. Genéralo con: samtools faidx {path}")
        self.idx = {}
        for linea in open(fai):
            c = linea.rstrip("\n").split("\t")
            if len(c) >= 5:
                self.idx[c[0]] = (int(c[1]), int(c[2]), int(c[3]), int(c[4]))
        self._fh = open(self.path, "rb")

    def __contains__(self, chrom):
        return chrom in self.idx

    def length(self, chrom):
        return self.idx[chrom][0]

    def fetch(self, chrom, start, end):
        """Devuelve la secuencia en mayúsculas. start y end son 0-based, semiabiertos.

        Se recorta a los límites del cromosoma; una petición entera fuera devuelve "".
        """
        if chrom not in self.idx:
            return ""
        n, off, bpl, bl = self.idx[chrom]
        start = max(0, start)
        end = min(n, end)
        if start >= end:
            return ""
        # byte donde empieza cada base: desplazamiento + líneas completas + resto
        def byte_de(pos):
            return off + (pos // bpl) * bl + (pos % bpl)
        self._fh.seek(byte_de(start))
        # se leen bytes de más por los saltos de línea y luego se limpian
        bruto = self._fh.read(byte_de(end) - byte_de(start))
        return bruto.replace(b"\n", b"").replace(b"\r", b"").upper().decode("ascii")

    def close(self):
        try:
            self._fh.close()
        except Exception:
            pass
