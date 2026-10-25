#!/usr/bin/env python3
"""Lee los cinco VCF del cromosoma 21 y produce alelos codificados 0/1 por individuo.

Recorre los cinco ficheros a la vez (todos están ordenados por posición) usando el VCF de
Mbuti como guía, porque es el único que trae el estado ancestral. Para cada posición de
Mbuti con CAnc declarado, busca esa misma posición en los cuatro arcaicos.

Codificación de salida, un array int8 por individuo:
    -1  sin dato        0  alelo ancestral        1  alelo derivado

"Derivado" significa distinto de CAnc, la base del ancestro humano-chimpancé del
alineamiento EPO de 6 primates (Ensembl v64), que viene en el propio VCF de Mbuti.

uso:
    python3 scripts/read_vcfs.py --dir ~/sima/public_vcfs/chr21 \
            --out results/chr21_alleles.npz [opciones]

Opciones principales:
    --mode pseudohaploid   un alelo al azar de los dos (por defecto; comparable con la simulación)
    --mode homozygous      solo genotipos homocigotos; los heterocigotos cuentan como sin dato
    --exclude-rm           descarta posiciones marcadas RM (repetición enmascarada)
    --min-dp N             exige profundidad mínima en los arcaicos (por defecto 0, sin filtro)
    --max-lines N          se detiene tras N líneas de Mbuti (para pruebas rápidas)
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np

SAMPLES = ["D5", "Vi3319", "D3", "Chag8"]      # arcaicos, en orden de salida
BASES = frozenset(b"ACGT")


def opener(path):
    """Devuelve un iterador de líneas en bytes, descomprimiendo con zcat (más rápido que gzip)."""
    p = subprocess.Popen(["zcat", str(path)], stdout=subprocess.PIPE, bufsize=1 << 20)
    return p, p.stdout


class ArchaicReader:
    """Avanza por un VCF arcaico de una sola muestra, en orden de posición."""

    __slots__ = ("proc", "fh", "pos", "ref", "alt", "gt", "dp", "eof", "name", "n_lines")

    def __init__(self, path, name):
        self.proc, self.fh = opener(path)
        self.name = name
        self.eof = False
        self.n_lines = 0
        self.pos = -1
        # _skip_header deja ya analizado el primer registro de datos: NO avanzar después,
        # o se perdería la primera posición del fichero.
        self._skip_header()

    def _skip_header(self):
        for line in self.fh:
            if not line.startswith(b"#"):
                self._parse(line)
                return
        self.eof = True
        self.pos = 1 << 62

    def _parse(self, line):
        f = line.rstrip(b"\n").split(b"\t")
        if len(f) < 10:
            self.eof = True
            return
        self.pos = int(f[1])
        self.ref = f[3]
        self.alt = f[4]
        s = f[9]
        c = s.find(b":")
        self.gt = s[:c] if c > 0 else s
        # DP es el segundo campo del formato GT:DP:...
        if c > 0:
            d = s.find(b":", c + 1)
            seg = s[c + 1:d] if d > 0 else s[c + 1:]
            self.dp = int(seg) if seg.isdigit() else 0
        else:
            self.dp = 0
        self.n_lines += 1

    def _next(self):
        line = self.fh.readline()
        if not line:
            self.eof = True
            self.pos = 1 << 62
            return
        self._parse(line)

    def seek(self, target):
        """Avanza hasta pos >= target. Devuelve True si quedó exactamente en target."""
        while not self.eof and self.pos < target:
            self._next()
        return (not self.eof) and self.pos == target

    def alleles(self):
        """Devuelve las dos bases del genotipo, o None si no es utilizable."""
        gt = self.gt
        if len(gt) != 3 or gt[1:2] not in (b"/", b"|"):
            return None
        a, b = gt[0:1], gt[2:3]
        if not (a.isdigit() and b.isdigit()):
            return None
        ia, ib = int(a), int(b)
        if ia > 1 or ib > 1:          # 1/2 y similares: multialélico, se descarta
            return None
        if ia or ib:
            if self.alt == b"." or b"," in self.alt or len(self.alt) != 1:
                return None
        if len(self.ref) != 1:
            return None
        alt = self.alt
        return (self.ref if ia == 0 else alt, self.ref if ib == 0 else alt)

    def close(self):
        try:
            self.proc.kill()
        except Exception:
            pass


def info_get(info, key):
    """Extrae key=valor de un campo INFO en bytes. Devuelve None si no está."""
    i = info.find(key)
    while i != -1:
        if i == 0 or info[i - 1:i] == b";":
            j = info.find(b";", i)
            return info[i + len(key):j if j != -1 else len(info)]
        i = info.find(key, i + 1)
    return None


def info_flag(info, flag):
    i = info.find(flag)
    while i != -1:
        before_ok = i == 0 or info[i - 1:i] == b";"
        after = info[i + len(flag):i + len(flag) + 1]
        if before_ok and after in (b"", b";"):
            return True
        i = info.find(flag, i + 1)
    return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--mode", choices=["pseudohaploid", "homozygous", "carrier"],
                    default="pseudohaploid")
    ap.add_argument("--seed", type=int, default=20261025)
    ap.add_argument("--exclude-rm", action="store_true")
    ap.add_argument("--min-dp", type=int, default=0)
    ap.add_argument("--max-lines", type=int, default=0)
    ap.add_argument("--progress", type=int, default=2_000_000)
    a = ap.parse_args()

    d = Path(a.dir).expanduser()
    for n in SAMPLES + ["Mbuti"]:
        if not (d / f"{n}_chr21.vcf.gz").exists():
            sys.exit(f"ERROR: falta {d}/{n}_chr21.vcf.gz")

    rng = np.random.default_rng(a.seed)
    arch = [ArchaicReader(d / f"{n}_chr21.vcf.gz", n) for n in SAMPLES]
    mproc, mfh = opener(d / "Mbuti_chr21.vcf.gz")

    pos_l, map20_l, rm_l = [], [], []
    code_l = {n: [] for n in SAMPLES + ["Mbuti"]}
    stats = dict(mbuti_lines=0, with_canc=0, canc_bad=0, rm_skipped=0,
                 mbuti_unusable=0, emitted=0, ref_mismatch=0)

    def code(pair, canc):
        """Codifica un genotipo como -1 (sin dato), 0 (ancestral) o 1 (derivado).

        pseudohaploid: un alelo al azar de los dos.
        homozygous:    solo homocigotos; los heterocigotos cuentan como sin dato.
        carrier:       derivado si PORTA al menos un alelo derivado (lectura literal de
                       "carries the derived allele" en Massilani et al. 2026).
        """
        if pair is None:
            return -1
        d0, d1 = pair[0] != canc, pair[1] != canc
        if a.mode == "carrier":
            return 1 if (d0 or d1) else 0
        if a.mode == "homozygous":
            return -1 if d0 != d1 else (1 if d0 else 0)
        return 1 if (d0 if rng.random() < 0.5 else d1) else 0

    for line in mfh:
        if line.startswith(b"#"):
            continue
        stats["mbuti_lines"] += 1
        if a.max_lines and stats["mbuti_lines"] > a.max_lines:
            break
        if a.progress and stats["mbuti_lines"] % a.progress == 0:
            print(f"  ...{stats['mbuti_lines']:,} líneas de Mbuti, "
                  f"{stats['emitted']:,} posiciones emitidas", file=sys.stderr, flush=True)

        f = line.rstrip(b"\n").split(b"\t")
        if len(f) < 10:
            continue
        info = f[7]
        canc = info_get(info, b"CAnc=")
        if canc is None:
            continue
        stats["with_canc"] += 1
        canc = canc.upper()
        if len(canc) != 1 or canc[0] not in BASES:
            stats["canc_bad"] += 1
            continue
        if a.exclude_rm and info_flag(info, b"RM"):
            stats["rm_skipped"] += 1
            continue

        pos = int(f[1])
        ref, alt = f[3], f[4]
        # genotipo de Mbuti
        s = f[9]
        c = s.find(b":")
        gt = s[:c] if c > 0 else s
        mb_pair = None
        if len(gt) == 3 and gt[1:2] in (b"/", b"|") and gt[0:1].isdigit() and gt[2:3].isdigit():
            ia, ib = int(gt[0:1]), int(gt[2:3])
            if ia <= 1 and ib <= 1 and len(ref) == 1 and (
                    (ia == 0 and ib == 0) or (alt != b"." and b"," not in alt and len(alt) == 1)):
                mb_pair = (ref if ia == 0 else alt, ref if ib == 0 else alt)
        mb = code(mb_pair, canc)
        if mb < 0:
            stats["mbuti_unusable"] += 1
            continue

        # los cuatro arcaicos en esta posición
        codes, ok_ref = [], True
        for r in arch:
            if not r.seek(pos):
                codes.append(-1)
                continue
            if r.ref != ref:
                ok_ref = False
            if a.min_dp and r.dp < a.min_dp:
                codes.append(-1)
                continue
            codes.append(code(r.alleles(), canc))
        if not ok_ref:
            stats["ref_mismatch"] += 1

        # anotaciones que se guardan para poder filtrar después sin releer los VCF.
        # Map20 solo aparece en una parte de las posiciones: se guarda -1 cuando falta,
        # de modo que el análisis pueda distinguir "ausente" de "valor bajo".
        m20 = info_get(info, b"Map20=")
        try:
            map20_l.append(float(m20) if m20 is not None else -1.0)
        except ValueError:
            map20_l.append(-1.0)
        rm_l.append(info_flag(info, b"RM"))

        pos_l.append(pos)
        code_l["Mbuti"].append(mb)
        for n, cd in zip(SAMPLES, codes):
            code_l[n].append(cd)
        stats["emitted"] += 1

    for r in arch:
        r.close()
    try:
        mproc.kill()
    except Exception:
        pass

    out = Path(a.out).expanduser()
    out.parent.mkdir(parents=True, exist_ok=True)
    arrays = {"pos": np.array(pos_l, dtype=np.int64),
              "map20": np.array(map20_l, dtype=np.float32),
              "rm": np.array(rm_l, dtype=bool)}
    for n in SAMPLES + ["Mbuti"]:
        arrays[n] = np.array(code_l[n], dtype=np.int8)
    arrays["meta_mode"] = np.array([a.mode])
    arrays["meta_seed"] = np.array([a.seed])
    arrays["meta_exclude_rm"] = np.array([a.exclude_rm])
    arrays["meta_min_dp"] = np.array([a.min_dp])
    np.savez_compressed(out, **arrays)

    print()
    print(f"Escrito: {out}")
    print(f"  modo {a.mode}, semilla {a.seed}, exclude_rm={a.exclude_rm}, min_dp={a.min_dp}")
    print(f"  líneas de Mbuti recorridas: {stats['mbuti_lines']:,}")
    print(f"  con CAnc declarado:         {stats['with_canc']:,}")
    print(f"  CAnc no utilizable:         {stats['canc_bad']:,}")
    if a.exclude_rm:
        print(f"  descartadas por RM:         {stats['rm_skipped']:,}")
    print(f"  genotipo de Mbuti no usable:{stats['mbuti_unusable']:,}")
    print(f"  POSICIONES EMITIDAS:        {stats['emitted']:,}")
    if stats["ref_mismatch"]:
        print(f"  AVISO: REF discrepante entre ficheros en {stats['ref_mismatch']:,} posiciones")
    print()
    m20 = arrays["map20"]
    print(f"  Map20 presente en {int((m20 >= 0).sum()):,} posiciones "
          f"({100.0*(m20 >= 0).mean():.1f}%); de ellas, con valor 1: "
          f"{int((m20 >= 1).sum()):,}")
    print(f"  marcadas RM: {int(arrays['rm'].sum()):,} "
          f"({100.0*arrays['rm'].mean():.1f}%)")
    print()
    print("  cobertura por individuo (posiciones con dato):")
    for n in SAMPLES + ["Mbuti"]:
        v = arrays[n]
        cov = int((v >= 0).sum())
        der = int((v == 1).sum())
        pct = 100.0 * cov / len(v) if len(v) else 0.0
        print(f"    {n:8s} {cov:10,} ({pct:5.1f}%)   derivados: {der:,}")


if __name__ == "__main__":
    main()
