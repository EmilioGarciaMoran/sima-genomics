#!/usr/bin/env python3
"""Paso 6: resume los resultados en results/REPORT.md (tablas leídas de los TSV, sin teclear números)."""
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def table(path, cols=None):
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    cols = cols or list(rows[0]) if rows else (cols or [])
    out = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    out += ["| " + " | ".join(r[c] for c in cols) + " |" for r in rows]
    return "\n".join(out), rows


def main(slices: Path, out: Path) -> None:
    lines = ["# Informe del pipeline nuclear reproducible", ""]
    man, _ = table(slices / "MANIFEST.tsv", ["file", "sample", "build", "reads", "md5"])
    lines += ["## Datos de entrada (slices/MANIFEST.tsv)", "", man, ""]

    t, rows = table(out / "direct_summary.tsv")
    lines += ["## Genotipado directo de las posiciones del catálogo", "",
              "Pseudohaploide (una lectura al azar, semilla fija). `reads_within_window`: posiciones con "
              "alguna lectura a ± WINDOW pb aunque ninguna cubra la base.", "", t, ""]
    called = []
    for r in csv.DictReader(open(out / "direct_genotypes.tsv"), delimiter="\t"):
        if r["state"] != "NA":
            flag = " (transición: posible daño)" if r["transition"] == "1" and r["state"] == "DER" else ""
            called.append(f"| {r['sample']} | {r['id']} | {r['gene']} | {r['chrom']}:{r['pos']} | "
                          f"{r['depth']} | {r['call']} | {r['state']}{flag} | {r['state_trimmed']} |")
    lines += ["### Posiciones con lectura", ""]
    if called:
        lines += ["| muestra | id | gen | posición | prof. | base | estado | estado sin extremos |",
                  "|---|---|---|---|---|---|---|---|"] + called
    else:
        lines += ["Ninguna posición del catálogo tiene lecturas encima en ninguna muestra."]
    lines += [""]

    t, _ = table(out / "gl_summary.tsv")
    lines += ["## Verosimilitudes genotípicas en las CDS", "",
              "Sitios con base no-REF: candidatos, no variantes; las transversiones lejos de los extremos "
              "de la lectura (última columna) son las menos afectadas por el daño.", "", t, ""]
    angsd = out / "angsd"
    if angsd.exists():
        lines += ["## ANGSD (paso 3b)", "", f"Versión: {(angsd / 'VERSION.txt').read_text().strip()}",
                  "", "Salidas en `results/angsd/` (`*.mafs.gz`, `*.glf.gz`).", ""]
    pwd = out / "pairwise"
    if (pwd / "distances.tsv").exists():
        lines += ["## Distancias por pares", "",
                  "Distancia p sobre el panel de transversiones (relativa al panel, no por pb). "
                  "EE por jackknife de bloques. La distancia de un espécimen de baja cobertura incluye su error: "
                  "compárala con la de un comparador al mismo simio.", ""]
        t, _ = table(pwd / "specimens_summary.tsv")
        lines += [t, ""]
        roles = {r["specimen"]: r["role"] for r in csv.DictReader(open(pwd / "specimens_summary.tsv"), delimiter="\t")}
        dist = list(csv.DictReader(open(pwd / "distances.tsv"), delimiter="\t"))
        for tgt in [n for n, r in roles.items() if r == "target"]:
            rs = sorted(((r["specimen_2"] if r["specimen_1"] == tgt else r["specimen_1"], r) for r in dist
                         if tgt in (r["specimen_1"], r["specimen_2"])), key=lambda x: float(x[1]["distance"]))
            lines += [f"### {tgt}", "", "| espécimen | sitios | distancia | EE |", "|---|---|---|---|"]
            lines += [f"| {o} | {r['sites_compared']} | {r['distance']} | {r['se_jackknife']} |" for o, r in rs]
            lines += [""]
        if (pwd / "d_stats.tsv").exists():
            t, _ = table(pwd / "d_stats.tsv")
            lines += ["### D(target; X, Y): alelos derivados compartidos", "",
                      "D > 0: el target comparte más alelos derivados con X que con Y. |Z| > 3: significativo.",
                      "", t, ""]
    (out / "REPORT.md").write_text("\n".join(lines))
    print(f"  informe: {out / 'REPORT.md'}")


if __name__ == "__main__":
    a = sys.argv[1:] + [None] * 2
    main(Path(a[0] or ROOT / "slices"), Path(a[1] or ROOT / "results"))
