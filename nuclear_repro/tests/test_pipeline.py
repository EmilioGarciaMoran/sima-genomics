"""Prueba de extremo a extremo: fixture sintético -> recorte -> run_all.sh -> resultados esperados.

    cd nuclear_repro && python3 -m pytest -q tests
"""
import csv
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pysam
import pytest

ROOT = Path(__file__).resolve().parent.parent
ENV = {**os.environ, "BUILD": "fixture"}


def sh(*cmd, check=True):
    return subprocess.run([str(c) for c in cmd], cwd=ROOT, env=ENV, capture_output=True, text=True, check=check)


def rows(path):
    return list(csv.DictReader(open(path), delimiter="\t"))


@pytest.fixture(scope="module")
def run(tmp_path_factory):
    d = tmp_path_factory.mktemp("fx")
    sh(sys.executable, "tests/make_fixture.py", d)
    full = d / "full"
    sh(sys.executable, "scripts/00_make_slices.py", "--ref", full / "ref.fa",
       "--bam", full / "SH_A.bam", full / "SH_B.bam", "--catalog", d / "catalog", "--out", d / "slices")
    r = sh("bash", "run_all.sh", d / "slices", d / "catalog", d / "results")
    return d, r


def test_slices_keep_only_target_regions(run):
    d, _ = run
    man = {r["file"]: r for r in rows(d / "slices" / "MANIFEST.tsv")}
    assert man["SH_A.bam"]["build"] == "fixture"
    assert int(man["SH_B.bam"]["reads"]) == 2           # las 3 lecturas de chr1:118000 quedan fuera
    with pysam.AlignmentFile(str(d / "slices" / "SH_B.bam")) as b:
        assert all(r.reference_start < 118000 for r in b)
    assert "chr1:48001-52001" in (d / "slices" / "regions.fa").read_text()


def test_direct_genotypes(run):
    d, _ = run
    g = {(r["sample"], r["id"]): r for r in rows(d / "results" / "direct_genotypes.tsv")}
    assert g["SH_A", "P1"]["state"] == "DER" and g["SH_A", "P1"]["depth"] == "10"
    assert g["SH_B", "P2"]["state"] == "ANC" and g["POOL", "P2"]["depth"] == "6"
    assert g["POOL", "P3"]["state"] == "NA" and int(g["POOL", "P3"]["reads_window"]) > 0
    assert g["POOL", "P4"]["state"] == "NA" and g["POOL", "P4"]["reads_window"] == "0"
    p5 = g["SH_A", "P5"]
    assert (p5["state"], p5["transition"], p5["state_trimmed"]) == ("DER", "1", "NA")
    s = {r["sample"]: r for r in rows(d / "results" / "direct_summary.tsv")}
    assert s["POOL"]["covered_ge1"] == "3" and s["POOL"]["DER_transversion"] == "1"


def test_genotype_likelihoods(run):
    d, _ = run
    sites = {(r["sample"], r["pos"]): r for r in rows(d / "results" / "gl_sites.tsv")}
    hom, het = sites["SH_A", "90101"], sites["SH_A", "100051"]
    assert hom["gt_bases"] == "CC" and hom["ref"] == "A"
    assert het["gt_bases"] == "GT" and het["best_gt"] == "Mm"
    assert float(het["GL_MM"]) < -3 and float(het["GL_mm"]) < -3   # het claramente mejor
    summ = {(r["sample"], r["gene"]): r for r in rows(d / "results" / "gl_summary.tsv")}
    assert summ["POOL", "GENE_EMPTY"]["bp_covered"] == "0"
    assert summ["SH_A", "GENE_HOM"]["sites_nonref_transversion_trimmed"] == "1"


def test_report_written(run):
    d, _ = run
    rep = (d / "results" / "REPORT.md").read_text()
    assert "Genotipado directo" in rep and "posible daño" in rep


def test_deterministic(run, tmp_path):
    d, _ = run
    sh("bash", "run_all.sh", d / "slices", d / "catalog", tmp_path)
    for f in ("direct_genotypes.tsv", "gl_sites.tsv", "gl_summary.tsv"):
        assert (tmp_path / f).read_bytes() == (d / "results" / f).read_bytes(), f


def test_build_mismatch_is_caught(run, tmp_path):
    """Un catálogo con REF cambiado (como al mezclar hg19/hg38) tiene que abortar el paso 1."""
    d, _ = run
    cat = tmp_path / "catalog"
    shutil.copytree(d / "catalog", cat)
    p = cat / "positions.tsv"
    p.write_text(p.read_text().replace("P1\tchr1\t50001\tC", "P1\tchr1\t50001\tG"))
    r = sh(sys.executable, "scripts/01_check_inputs.py", d / "slices", cat, check=False)
    assert r.returncode != 0 and "REF distinto" in r.stderr


def test_tampered_slice_is_caught(run, tmp_path):
    d, _ = run
    sl = tmp_path / "slices"
    shutil.copytree(d / "slices", sl)
    with open(sl / "regions.fa", "a") as fh:
        fh.write("\n")
    r = sh(sys.executable, "scripts/01_check_inputs.py", sl, d / "catalog", check=False)
    assert r.returncode != 0 and "md5" in r.stderr
