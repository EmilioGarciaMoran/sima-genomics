#!/usr/bin/env python3
"""Prueba de respuesta conocida para el cálculo de D de analyze_sima.py.

Construye sitios y alelos donde la D correcta se conoce por aritmética.
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent.parent / "branch_counts" / "scripts"))
import importlib.util
spec = importlib.util.spec_from_file_location("an", HERE / "analyze_sima.py")
an = importlib.util.module_from_spec(spec)
sys.argv = ["x", "--sitios", "/dev/null", "--alelos", "/dev/null"]
try:
    spec.loader.exec_module(an)
except SystemExit:
    pass

fails = []
def check(name, got, want, tol=1e-9):
    ok = abs(got - want) <= tol if isinstance(want, float) else got == want
    print(f"  {'ok  ' if ok else 'FALLA'} {name}: obtenido {got}, esperado {want}")
    if not ok: fails.append(name)

print("1. D con recuentos conocidos")
# 60 sitios ABBA (Sima con Vindija) y 20 BABA -> D = (60-20)/80 = 0.5
n = 80
d3 = np.array([0]*60 + [1]*20, dtype=np.int8)
vi = np.array([1]*60 + [0]*20, dtype=np.int8)
sima = np.ones(n, dtype=np.int8)
mb = np.zeros(n, dtype=np.int8)
blocks = np.arange(n) % 40          # 40 bloques
(D, ee, nn), ab, ba = an.calcula_D(sima, d3, vi, mb, blocks)
check("ABBA", ab, 60); check("BABA", ba, 20); check("n", nn, 80)
check("D", round(D, 6), 0.5)
check("el EE es finito y positivo", bool(np.isfinite(ee) and ee > 0), True)

print("\n2. D = 0 con recuentos iguales")
d3 = np.array([0]*40 + [1]*40, dtype=np.int8)
vi = np.array([1]*40 + [0]*40, dtype=np.int8)
(D, ee, nn), ab, ba = an.calcula_D(sima, d3, vi, mb, blocks)
check("D nula", round(D, 6), 0.0); check("n", nn, 80)

print("\n3. D negativa si Sima va con Denisova")
d3 = np.array([0]*20 + [1]*60, dtype=np.int8)
vi = np.array([1]*20 + [0]*60, dtype=np.int8)
(D, _, _), ab, ba = an.calcula_D(sima, d3, vi, mb, blocks)
check("D = -0.5", round(D, 6), -0.5)

print("\n4. Solo cuentan los sitios donde Sima lleva el derivado")
d3 = np.array([0]*60 + [1]*20, dtype=np.int8)
vi = np.array([1]*60 + [0]*20, dtype=np.int8)
s2 = sima.copy(); s2[:30] = 0        # 30 ABBA pasan a ancestral
(D, _, nn), ab, ba = an.calcula_D(s2, d3, vi, mb, blocks)
check("ABBA baja a 30", ab, 30); check("n baja a 50", nn, 50)
check("D = (30-20)/50", round(D, 6), 0.2)

print("\n5. Sin dato en Sima no cuenta")
s3 = sima.copy(); s3[:30] = -1
(D, _, nn), ab, _ = an.calcula_D(s3, d3, vi, mb, blocks)
check("ABBA baja a 30", ab, 30); check("n = 50", nn, 50)

print("\n6. Mbuti derivado excluye el sitio")
mb2 = mb.copy(); mb2[:60] = 1
(D, _, nn), ab, ba = an.calcula_D(sima, d3, vi, mb2, blocks)
check("los 60 ABBA quedan fuera", ab, 0); check("quedan los 20 BABA", ba, 20)

print("\n7. Sitios no informativos (D3 y Vi iguales) no cuentan")
d3b = np.ones(n, dtype=np.int8); vib = np.ones(n, dtype=np.int8)
(D, _, nn), ab, ba = an.calcula_D(sima, d3b, vib, mb, blocks)
check("n = 0", nn, 0)

print("\n8. El jackknife se acerca al binomial con bloques homogéneos")
rng = np.random.default_rng(1); N = 4000
d3 = np.zeros(N, np.int8); vi = np.zeros(N, np.int8)
es_abba = rng.random(N) < 0.75
d3[~es_abba] = 1; vi[es_abba] = 1
sima = np.ones(N, np.int8); mb = np.zeros(N, np.int8)
blocks = np.arange(N) % 200
(D, ee, nn), ab, ba = an.calcula_D(sima, d3, vi, mb, blocks)
eb = np.sqrt((1 - D**2) / nn)
check("cociente EE jackknife / binomial", round(ee/eb, 2), 1.0, 0.25)

print()
if fails:
    print(f"RESULTADO: {len(fails)} prueba(s) fallidas: {', '.join(fails)}")
    sys.exit(1)
print("RESULTADO: todas las pruebas pasan")
