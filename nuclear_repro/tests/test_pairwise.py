"""Distancias por pares sobre un árbol simulado: ventanas -> nube (panel + vectores por ventana)
-> matriz -> estadísticos. Incluye la acumulación: ampliar ventanas o especímenes solo calcula lo que falta.

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


def windows(d, *chroms):
    """Ventanas de 10 kb consecutivas (en datos reales serían p. ej. 100.000-200.000 de cada cromosoma)."""
    return sh(sys.executable, "cloud/make_windows.py", "--fai", d / "data" / "ref.fa.fai", "--chroms", *chroms,
              "--start", 0, "--size", 10000, "--every", 10000, "--out", d / "pairwise" / "windows.bed")


def cloud(d, check=True, **env):
    return sh("bash", "cloud/run_vectors.sh", d / "data", d / "data" / "ref.fa", 2, check=check,
              env={"PAIRWISE_DIR": str(d / "pairwise"), "HEARTBEAT_MIN": "1", **env})


@pytest.fixture(scope="module")
def pw(tmp_path_factory):
    d = tmp_path_factory.mktemp("pw")
    sh(sys.executable, "tests/make_fixture_pairwise.py", d)
    windows(d, "chr1", "chr2")
    cloud(d)
    sh(sys.executable, "scripts/05_pairwise.py", d / "pairwise", d / "results")
    return d


def dist(d):
    return {frozenset((r["specimen_1"], r["specimen_2"])): r for r in rows(d / "results" / "distances.tsv")}


def dstats(d):
    return {(r["target"], r["X"], r["Y"]): r for r in rows(d / "results" / "d_stats.tsv")}


def test_vectors_built_for_every_specimen(pw):
    wins = sorted((pw / "pairwise" / "windows").glob("chr*"))
    assert len(wins) == 24                                       # 12 ventanas de 10 kb por cromosoma
    for w in wins:
        vecs = sorted(p.name for p in (w / "vectors").glob("*.vec.gz"))
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
    w = sorted((pw / "pairwise" / "windows").glob("chr*"))[3]
    for name, flag, src in (("Sima", "--bam", "Sima.bam"), ("Nea", "--vcf", "Nea.vcf.gz")):
        sh(sys.executable, "cloud/10_specimen_vector.py", "--panel", w / "panel.tsv.gz",
           "--name", name, flag, pw / "data" / src, "--out", tmp_path)
        assert (tmp_path / f"{name}.vec.gz").read_bytes() == (w / "vectors" / f"{name}.vec.gz").read_bytes()


def snapshot(d):
    return {p.relative_to(d): p.read_bytes() for p in (d / "pairwise" / "windows").rglob("*.gz")}


def test_accumulation_matches_single_run(pw, tmp_path):
    """Primero chr1, luego se añade chr2: lo ya hecho no se toca y el resultado es idéntico al de una vez."""
    d = tmp_path
    shutil.copytree(pw / "data", d / "data")
    shutil.copytree(pw / "pairwise", d / "pairwise", ignore=shutil.ignore_patterns("windows*", "*.log"))
    windows(d, "chr1")
    cloud(d)
    first = snapshot(d)
    assert first and all("chr1_" in str(p) for p in first)
    assert "12 ventanas nuevas, 24 en total" in windows(d, "chr2").stdout
    r = cloud(d)
    assert "trabajos: 12 paneles + 84 vectores" in r.stdout    # solo las ventanas nuevas
    second = snapshot(d)
    assert all(second[p] == b for p, b in first.items())      # chr1 intacto, byte a byte
    sh(sys.executable, "scripts/05_pairwise.py", d / "pairwise", d / "results")
    for f in ("distances.tsv", "d_stats.tsv", "derived_sharing.tsv", "windows_summary.tsv"):
        assert (d / "results" / f).read_bytes() == (pw / "results" / f).read_bytes(), f


def test_new_or_missing_vector_is_computed_alone(pw, tmp_path):
    d = tmp_path
    shutil.copytree(pw / "data", d / "data")
    shutil.copytree(pw / "pairwise", d / "pairwise", ignore=shutil.ignore_patterns("*.log"))
    target = sorted((d / "pairwise" / "windows").glob("chr*"))[5] / "vectors" / "Gorilla.vec.gz"
    target.unlink()
    r = cloud(d)
    assert "trabajos: 0 paneles + 1 vectores" in r.stdout
    assert target.exists()


def test_pending_windows_are_ignored_and_listed(pw, tmp_path):
    d = tmp_path
    shutil.copytree(pw / "pairwise", d / "pairwise", ignore=shutil.ignore_patterns("*.log"))
    with open(d / "pairwise" / "windows.bed", "a") as fh:
        fh.write("chr2\t120000\t130000\tchr2_000120001_000130000\n")
    r = sh(sys.executable, "scripts/05_pairwise.py", d / "pairwise", d / "res")
    assert "1 ventanas de windows.bed sin procesar" in r.stdout
    assert (d / "res" / "distances.tsv").read_bytes() == (pw / "results" / "distances.tsv").read_bytes()


def test_vector_from_other_panel_is_rejected(pw, tmp_path):
    p = tmp_path / "pairwise"
    shutil.copytree(pw / "pairwise", p)
    panel = sorted((p / "windows").glob("chr*"))[0] / "panel.tsv.gz"
    with gzip.open(panel, "rt") as fh:
        lines = fh.read().splitlines()
    with gzip.open(panel, "wt") as fh:
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
