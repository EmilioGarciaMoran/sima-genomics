#!/usr/bin/env python3
"""Pruebas de respuesta conocida para extract_alleles.py.

Construye un conjunto de sitios y un SAM diminutos, con la respuesta conocida por
construcción.

    python3 scripts/test_extract_alleles.py
"""
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
fails = []


def check(name, got, want):
    ok = got == want
    print(f"  {'ok  ' if ok else 'FALLA'} {name}: obtenido {got}, esperado {want}")
    if not ok:
        fails.append(name)


# Sitios: cromosoma 1, posiciones 100..600 de cien en cien.
# En todos, ancestral = A, derivado = G.
POS = np.array([100, 200, 300, 400, 500, 600], dtype=np.int64)
SITIOS = dict(
    pos=POS,
    chrom=np.ones(len(POS), dtype=np.int8),
    canc=np.frombuffer(b"AAAAAA", dtype=np.uint8).copy(),
    der=np.frombuffer(b"GGGGGG", dtype=np.uint8).copy(),
)
IDX = {int(p): i for i, p in enumerate(POS)}


def rec(nombre, pos, cigar, seq, cr="1", flag=0):
    return f"{nombre}\t{flag}\t{cr}\t{pos}\t60\t{cigar}\t*\t0\t0\t{seq}\t{'I'*len(seq)}\n"


def run(d, sam_texto, *args):
    sam = d / "t.sam"
    sam.write_text(sam_texto)
    out = d / "o.npz"
    r = subprocess.run([sys.executable, str(HERE / "extract_alleles.py"),
                        "--sitios", str(d / "s.npz"), "--bam", str(sam), "--sam",
                        "--out", str(out), *args], capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stderr, file=sys.stderr)
        raise SystemExit("extract_alleles.py terminó con error")
    z = np.load(out, allow_pickle=False)
    return z["codigo"], z["cobertura"], r.stderr


with tempfile.TemporaryDirectory() as td:
    d = Path(td)
    np.savez_compressed(d / "s.npz", **SITIOS)

    print("1. Codificación básica")
    # lectura de 10 pb desde la posicion 100: cubre 100..109, solo el sitio 100
    sam = (rec("ancestral", 100, "10M", "A" + "C" * 9)
           + rec("derivado", 200, "10M", "G" + "C" * 9)
           + rec("tercero", 300, "10M", "T" + "C" * 9)
           + rec("enmascarado", 400, "10M", "N" + "C" * 9))
    cod, cob, err = run(d, sam)
    check("sitio 100 ancestral", int(cod[IDX[100]]), 0)
    check("sitio 200 derivado", int(cod[IDX[200]]), 1)
    check("sitio 300 tercer alelo -> sin dato", int(cod[IDX[300]]), -1)
    check("sitio 400 enmascarado -> sin dato", int(cod[IDX[400]]), -1)
    check("sitio 500 sin cubrir", int(cod[IDX[500]]), -1)
    check("cobertura del sitio 100", int(cob[IDX[100]]), 1)
    check("cobertura del sitio 400 (N no cuenta)", int(cob[IDX[400]]), 0)

    print("\n2. Una lectura puede cubrir varios sitios")
    sam = rec("larga", 100, "501M", "A" + "C" * 99 + "G" + "C" * 99 + "A"
              + "C" * 99 + "G" + "C" * 99 + "A" + "C" * 100)
    cod, _, _ = run(d, sam)
    check("cinco sitios codificados", [int(cod[IDX[p]]) for p in (100, 200, 300, 400, 500)],
          [0, 1, 0, 1, 0])
    check("el sexto queda fuera", int(cod[IDX[600]]), -1)

    print("\n3. Deleción: el sitio tapado queda sin dato")
    # 5M 10D 486M desde 100: cubre 100-104, salta 105-114, sigue en 115
    sam = rec("del", 100, "5M10D486M", "A" + "C" * 4 + "C" * 85 + "G" + "C" * 400)
    cod, _, _ = run(d, sam)
    check("sitio 100 sigue codificado", int(cod[IDX[100]]), 0)
    check("sitio 200 codificado tras la delecion", int(cod[IDX[200]]), 1)

    print("\n4. Inserción: no desplaza las posiciones de referencia")
    # 50M 5I 451M: la insercion consume secuencia pero no referencia
    # La secuencia se construye a partir del CIGAR para no equivocarse con los indices:
    # 50M (ref 100-149, seq 0-49), 5I (seq 50-54), 451M (ref 150-600, seq 55-505).
    def construir(cigar, ini, bases):
        """bases: {posicion_de_referencia: base}. El resto, C; las inserciones, T."""
        import re as _re
        s, r, q_ = [], ini, 0
        for num, op in _re.findall(r"(\d+)([MIDNSHP=X])", cigar):
            num = int(num)
            if op in ("M", "=", "X"):
                for i in range(num):
                    s.append(bases.get(r + i, "C"))
                r += num; q_ += num
            elif op in ("I", "S"):
                s.append("T" * num); q_ += num
            elif op in ("D", "N"):
                r += num
        return "".join(s)

    seq = construir("50M5I451M", 100, {100: "A", 200: "G", 300: "A"})
    assert len(seq) == 506, len(seq)
    sam = rec("ins", 100, "50M5I451M", seq)
    cod, _, _ = run(d, sam)
    check("sitio 100 correcto", int(cod[IDX[100]]), 0)
    check("sitio 200 correcto pese a la insercion", int(cod[IDX[200]]), 1)
    check("sitio 300 tambien", int(cod[IDX[300]]), 0)

    print("\n5. Recorte blando: tampoco desplaza")
    sam = rec("rec", 100, "5S10M", "TTTTT" + "A" + "C" * 9)
    cod, _, _ = run(d, sam)
    check("sitio 100 correcto con 5S delante", int(cod[IDX[100]]), 0)

    print("\n6. Muestreo pseudohaploide con varias lecturas")
    sam = "".join(rec(f"r{i}", 100, "10M", "A" + "C" * 9) for i in range(5)) \
        + "".join(rec(f"s{i}", 100, "10M", "G" + "C" * 9) for i in range(5))
    cod, cob, _ = run(d, sam)
    check("cobertura 10", int(cob[IDX[100]]), 10)
    check("el codigo es 0 o 1, no un consenso", int(cod[IDX[100]]) in (0, 1), True)
    c1, _, _ = run(d, sam, "--seed", "7")
    c2, _, _ = run(d, sam, "--seed", "7")
    check("misma semilla, mismo resultado", int(c1[IDX[100]]) == int(c2[IDX[100]]), True)
    # con semillas distintas sobre 50/50 debe salir cada valor alguna vez
    vistos = {int(run(d, sam, "--seed", str(k))[0][IDX[100]]) for k in range(1, 30)}
    check("con semillas distintas salen ambos alelos", sorted(vistos), [0, 1])

    print("\n6 bis. Un CIGAR incoherente se descarta sin tumbar el análisis")
    sam = (rec("mala", 100, "50M", "A" + "C" * 9)        # el CIGAR dice 50, la seq tiene 10
           + rec("buena", 200, "10M", "G" + "C" * 9))
    cod, _, err = run(d, sam)
    check("la mala no codifica nada", int(cod[IDX[100]]), -1)
    check("la buena si", int(cod[IDX[200]]), 1)
    check("se avisa del descarte", "CIGAR incoherente" in err, True)

    print("\n7. Otro cromosoma no contamina")
    sam = rec("otro_cr", 100, "10M", "G" + "C" * 9, cr="2")
    cod, _, _ = run(d, sam)
    check("sitio 100 del cr1 sigue sin dato", int(cod[IDX[100]]), -1)

    print("\n8. Cromosomas no numéricos se ignoran sin fallar")
    sam = rec("x", 100, "10M", "G" + "C" * 9, cr="X") + rec("mt", 100, "10M", "G" + "C" * 9, cr="MT")
    cod, _, _ = run(d, sam)
    check("ningun sitio codificado", int((cod >= 0).sum()), 0)

    print("\n9. Las estadísticas cuadran")
    sam = (rec("a", 100, "10M", "A" + "C" * 9)
           + rec("b", 200, "10M", "G" + "C" * 9)
           + rec("c", 300, "10M", "N" + "C" * 9))
    cod, _, err = run(d, sam)
    def val(etq):
        for l in err.splitlines():
            if etq in l:
                return int(l.split(":")[1].strip().replace(",", "").split()[0])
        return None
    check("lecturas", val("lecturas procesadas"), 3)
    check("bases alineadas", val("bases alineadas"), 30)
    check("caen en un sitio", val("caen en un sitio"), 3)
    check("enmascaradas", val("de ellas enmascaradas"), 1)
    check("sitios con dato", val("SITIOS CON DATO"), 2)
    check("derivados", val("derivados"), 1)

print()
if fails:
    print(f"RESULTADO: {len(fails)} prueba(s) fallidas: {', '.join(fails)}")
    sys.exit(1)
print("RESULTADO: todas las pruebas pasan")
