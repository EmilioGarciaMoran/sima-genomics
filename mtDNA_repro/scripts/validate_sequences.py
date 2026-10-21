#!/usr/bin/env python3
"""Paso 0 (validación): comprueba que cada secuencia de sequences/ es lo que dice ser.

Motivo: en la fase 2 del proyecto, 10 accesiones escritas a mano correspondían a hongos,
una polilla, un oomiceto, un parásito de peces, un cloroplasto de alga y un virus, y nadie
lo comprobó antes de alinear (ver CORRECTION.md). Este script aborta si algo así se repite.

Comprobaciones por taxón de accessions.tsv:
  1. existe sequences/<taxon>.fasta con la accesión esperada en la cabecera
  2. si viene de GenBank: el título original contiene 'mitochondri', 'complete' y el
     organismo esperado (columna expected_title_regex)
  3. longitud de genoma mitocondrial completo (15.000-17.500 pb)
  4. composición: >= 95 % de A/C/G/T/N y sin caracteres extraños
  5. el hash de la secuencia coincide con sequences/GENBANK_VERIFIED.tsv (si existe)

uso: validate_sequences.py accessions.tsv sequences/   -> código de salida 0 si todo es correcto
"""
import hashlib
import re
import sys
from pathlib import Path


def main(acc_tsv, seq_dir):
    seq_dir = Path(seq_dir)
    rows = [l.split("\t") for l in Path(acc_tsv).read_text().splitlines()[1:] if l.strip()]
    verified = {}
    vf = seq_dir / "GENBANK_VERIFIED.tsv"
    if vf.exists():
        for l in vf.read_text().splitlines()[1:]:
            if l.strip():
                t, a, h = l.split("\t")[:3]
                verified[t] = (a, h)

    errors = []
    for r in rows:
        taxon, acc, regex = r[0], r[1], r[5] if len(r) > 5 else ""
        f = seq_dir / f"{taxon}.fasta"
        if not f.exists():
            errors.append(f"{taxon}: falta {f}")
            continue
        lines = f.read_text().splitlines()
        header, seq = lines[0], "".join(lines[1:]).upper()

        if acc not in header:
            errors.append(f"{taxon}: la cabecera no contiene la accesión {acc}: {header[:80]}")
        if "source=GenBank" in header:
            title = header.split("title=", 1)[1] if "title=" in header else ""
            low = title.lower()
            if "mitochondri" not in low or "complete" not in low:
                errors.append(f"{taxon}: el título GenBank no es un genoma mitocondrial completo: {title[:90]}")
            if regex and not re.search(regex, title):
                errors.append(f"{taxon}: el título GenBank no corresponde al organismo esperado ({regex}): {title[:90]}")
        if not 15000 <= len(seq) <= 17500:
            errors.append(f"{taxon}: longitud {len(seq)} pb fuera de rango para un mtDNA completo")
        acgtn = sum(seq.count(c) for c in "ACGTN")
        if len(seq) == 0 or acgtn / len(seq) < 0.95:
            errors.append(f"{taxon}: composición anómala ({acgtn}/{len(seq)} A/C/G/T/N)")
        bad = set(seq) - set("ACGTNRYKMSWBDHV")
        if bad:
            errors.append(f"{taxon}: caracteres no válidos {sorted(bad)}")
        h = hashlib.sha256(seq.encode()).hexdigest()[:16]
        if taxon in verified and verified[taxon] != (acc, h):
            errors.append(f"{taxon}: hash {h} distinto del verificado {verified[taxon][1]}")
        status = "OK" if not any(e.startswith(taxon + ":") for e in errors) else "FALLO"
        print(f"{status:5} {taxon:11} {acc:12} {len(seq):6} pb  hash={h}")

    if errors:
        print("\nVALIDACIÓN FALLIDA:", *errors, sep="\n  - ", file=sys.stderr)
        sys.exit(1)
    print("Validación correcta: todas las secuencias son genomas mitocondriales completos de los taxones esperados.")


if __name__ == "__main__":
    main(*sys.argv[1:3])
