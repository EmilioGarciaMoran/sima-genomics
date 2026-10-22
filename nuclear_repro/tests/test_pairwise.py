"""Distancias por pares sobre un árbol simulado: nube (panel + vectores) -> matriz -> estadísticos.

    cd nuclear_repro && python3 -m pytest -q tests/test_pairwise.py
"""
import csv
import gzip
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent


def sh(*cmd, env=None, check=True):
    return subprocess.run([str(c) for c in cmd], cwd=ROOT, env={**os.environ, **(env or {})},
                          capture_output=True, text=True, check=check)


def rows(path):
    return list(csv.DictReader(open(path), delimiter="\t"))


@pytest.fixture(scope="module")
def pw(tmp_path_factory):
    d = tmp_path_factory.mktemp("pw")
    sh(sys.executable, "tests/make_fixture_pairwise.py", d)
    sh("bash", "cloud/run_vectors.sh", d / "data", d / "data" / "ref.fa", 2,
       env={"PAIRWISE_DIR": str(d / "pairwise"), "CHROMS": "chr1 chr2"})
    sh(sys.executable, "scripts/05_pairwise.py", d / "pairwise", d / "results", env={"BLOCK_SIZE": "10000"})
    return d


def dist(d):
    return {frozenset((r["specimen_1"], r["specimen_2"])): r for r in rows(d / "results" / "distances.tsv")}


def dstats(d):
    return {(r["target"], r["X"], r["Y"]): r for r in rows(d / "results" / "d_stats.tsv")}


def test_vectors_built_for_every_specimen(pw):
    vecs = sorted(p.name for p in (pw / "pairwise" / "vectors").glob("*.vec.gz"))
    assert vecs == sorted(f"{n}.vec.gz" for n in ("Sima", "Sima_true", "Nea", "Den", "Human", "Chimp", "Gorilla"))
    summ = {r["specimen"]: float(r["frac_called"]) for r in rows(pw / "results" / "specimens_summary.tsv")}
    assert 0.5 < summ["Sima"] < 0.85            # ~1x: alrededor de 1 - e^-1 de los sitios
    assert summ["Nea"] == 1.0


def test_distances_recover_the_tree(pw):
    D = dist(pw)
    d = lambda a, b: float(D[frozenset((a, b))]["distance"])  # noqa: E731
    assert d("Sima", "Nea") < d("Sima", "Den") < d("Sima", "Human") < d("Sima", "Chimp")
    assert d("Nea", "Den") < d("Nea", "Human")
    assert 0 < d("Sima", "Sima_true") < 0.02     # el error de secuenciación, visible y pequeño


def test_d_statistics(pw):
    S = dstats(pw)
    nd = S["Sima", "Nea", "Den"]
    assert float(nd["D"]) > 0 and float(nd["Z"]) > 3          # Sima más cerca de Nea que de Den
    assert float(S["Sima", "Den", "Human"]["Z"]) > 3           # y más cerca de Den que del humano
    assert not any("Gorilla" in k for k in S)                  # los simios no entran en D


def test_low_coverage_matches_truth(pw):
    """Con ~1x y 2% de error, D y los alelos derivados compartidos de Sima coinciden con los verdaderos."""
    S = dstats(pw)
    a, b = S["Sima", "Nea", "Den"], S["Sima_true", "Nea", "Den"]
    tol = 3 * math.hypot(float(a["se_jackknife"]), float(b["se_jackknife"]))
    assert abs(float(a["D"]) - float(b["D"])) < tol
    ds = {(r["target"], r["comparator"]): r for r in rows(pw / "results" / "derived_sharing.tsv")}
    for x in ("Nea", "Den", "Human"):
        s, t = ds["Sima", x], ds["Sima_true", x]
        tol = 3 * math.hypot(float(s["se_jackknife"]), float(t["se_jackknife"]))
        assert abs(float(s["fraction"]) - float(t["fraction"])) < tol, x


def test_vectors_are_deterministic(pw, tmp_path):
    for name, flag, src in (("Sima", "--bam", "Sima.bam"), ("Nea", "--vcf", "Nea.vcf.gz")):
        sh(sys.executable, "cloud/10_specimen_vector.py", "--panel", pw / "pairwise" / "panel.tsv.gz",
           "--name", name, flag, pw / "data" / src, "--out", tmp_path)
        assert (tmp_path / f"{name}.vec.gz").read_bytes() == (pw / "pairwise" / "vectors" / f"{name}.vec.gz").read_bytes()


def test_vector_from_other_panel_is_rejected(pw, tmp_path):
    p = tmp_path / "pairwise"
    shutil.copytree(pw / "pairwise", p)
    with gzip.open(p / "panel.tsv.gz", "rt") as fh:
        lines = fh.read().splitlines()
    with gzip.open(p / "panel.tsv.gz", "wt") as fh:
        fh.write("\n".join(lines[:-1]) + "\n")                  # un sitio menos
    r = sh(sys.executable, "scripts/05_pairwise.py", p, tmp_path / "res", check=False)
    assert r.returncode != 0 and "otro panel" in r.stderr


def test_allele_encoding():
    """Un alelo que no es REF ni ALT (error, daño, tercer alelo) es N, nunca REF: si no, los errores de
    un espécimen de baja cobertura lo acercarían a la referencia (al humano). Este sesgo (~1% de los
    sitios en el caso simulado) cae dentro del ruido de las pruebas de extremo a extremo, por eso va aquí."""
    sys.path.insert(0, str(ROOT / "scripts"))
    from pairwise_io import encode_allele
    assert [encode_allele(b, "A", "C") for b in ("A", "C", "G", "T", None)] == ["0", "1", "N", "N", "N"]
