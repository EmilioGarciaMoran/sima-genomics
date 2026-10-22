"""Modos de fallo de la ejecución en la nube (cloud/run_vectors.sh): que nada se pierda ni se dispare.

Se usa un 'aws' falso en el PATH que solo registra las llamadas, para comprobar que la subida a S3
ocurre también cuando hay fallos, hora límite o una interrupción.
"""
import os
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent


@pytest.fixture(scope="module")
def base(tmp_path_factory):
    d = tmp_path_factory.mktemp("cloud")
    subprocess.run([sys.executable, "tests/make_fixture_pairwise.py", str(d)], cwd=ROOT, check=True,
                   capture_output=True)
    subprocess.run([sys.executable, "cloud/make_windows.py", "--fai", str(d / "data" / "ref.fa.fai"),
                    "--chroms", "chr1", "chr2", "--start", "0", "--size", "10000", "--every", "10000",
                    "--out", str(d / "pairwise" / "windows.bed")], cwd=ROOT, check=True, capture_output=True)
    return d


@pytest.fixture
def run(base, tmp_path):
    shutil.copytree(base / "data", tmp_path / "data")
    shutil.copytree(base / "pairwise", tmp_path / "pairwise")
    (tmp_path / "bin").mkdir()
    aws = tmp_path / "bin" / "aws"
    aws.write_text(f'#!/bin/sh\necho "aws $*" >> {tmp_path}/aws.log\n')
    aws.chmod(0o755)

    def go(*, check=False, popen=False, **env):
        e = {**os.environ, "PAIRWISE_DIR": str(tmp_path / "pairwise"), "BUCKET": "s3://falso",
             "PATH": f"{tmp_path / 'bin'}:{os.environ['PATH']}", **env}
        cmd = ["bash", "cloud/run_vectors.sh", str(tmp_path / "data"), str(tmp_path / "data" / "ref.fa"), "2"]
        if popen:
            return subprocess.Popen(cmd, cwd=ROOT, env=e, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                    text=True, start_new_session=True)
        return subprocess.run(cmd, cwd=ROOT, env=e, capture_output=True, text=True, check=check, timeout=240)

    go.d = tmp_path
    go.aws = lambda: (tmp_path / "aws.log").read_text() if (tmp_path / "aws.log").exists() else ""
    go.vectors = lambda: sorted((tmp_path / "pairwise" / "windows").rglob("*.vec.gz"))
    return go


def test_plan_runs_nothing(run):
    r = run(PLAN="1", check=True)
    assert "trabajos: 24 paneles + 168 vectores" in r.stdout and "no se ejecuta nada" in r.stdout
    assert run.vectors() == [] and run.aws() == ""                 # ni cálculo ni S3


def test_pilot_processes_only_max_windows_and_uploads(run):
    r = run(MAX_WINDOWS="1", check=True)
    assert "en esta ejecución: 1 (MAX_WINDOWS=1)" in r.stdout
    assert len(run.vectors()) == 7
    log = run.aws()
    assert ".probe" in log.splitlines()[0]                          # S3 comprobado antes de calcular
    assert "sync" in log and "windows/" in log


def test_one_bad_input_does_not_stop_the_rest(run):
    (run.d / "data" / "Gorilla.fa").write_text(">chr1\nACGT\n")      # FASTA corrupto (el .fai no cuadra)
    r = run(MAX_WINDOWS="3")
    assert r.returncode == 1                                         # se avisa del fallo...
    assert "FALLIDOS" in r.stdout and r.stdout.count("/Gorilla") >= 3
    assert len(run.vectors()) == 3 * 6                               # ...pero el resto está hecho
    assert "sync" in run.aws()                                       # y subido


def test_deadline_leaves_jobs_pending(run):
    r = run(MAX_WINDOWS="1", MAX_HOURS="0", check=True)
    assert "sin empezar por la hora límite: 8" in r.stdout and run.vectors() == []
    r = run(MAX_WINDOWS="1", check=True)                             # la siguiente ejecución lo recoge
    assert len(run.vectors()) == 7


def test_missing_input_fails_before_any_work(run):
    (run.d / "data" / "Den.vcf.gz").unlink()
    r = run()
    assert r.returncode == 1 and "faltan ficheros de entrada" in r.stderr
    assert run.vectors() == [] and run.aws() == ""


def test_interrupt_still_uploads_and_resume_finishes(run):
    p = run(popen=True)
    deadline = time.time() + 60
    while not list((run.d / "pairwise" / "windows").glob("*/panel.tsv.gz")) and time.time() < deadline:
        time.sleep(0.1)
    os.killpg(p.pid, signal.SIGINT)                                  # Ctrl-C a todo el grupo, como un terminal
    time.sleep(0.3)
    os.killpg(p.pid, signal.SIGINT)                                  # y un segundo, impaciente
    out, _ = p.communicate(timeout=120)
    assert p.returncode == 130 and "interrumpido" in out and "Ejecución" in out
    assert "sync" in run.aws()
    r = run(check=True)
    assert "0 fallidos" in r.stdout and len(run.vectors()) == 24 * 7


def test_second_concurrent_run_is_refused(run):
    p = run(popen=True)
    deadline = time.time() + 60
    while not (run.d / "pairwise" / "runs").exists() and time.time() < deadline:   # ya tiene el cerrojo
        time.sleep(0.05)
    r = run()
    assert r.returncode == 1 and "ya hay otra ejecución en marcha" in r.stderr
    p.communicate(timeout=240)
    assert p.returncode == 0
