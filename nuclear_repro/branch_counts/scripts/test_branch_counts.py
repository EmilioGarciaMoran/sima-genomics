#!/usr/bin/env python3
"""Pruebas de respuesta conocida para branch_counts.py.

Cada caso se construye a mano, de modo que la respuesta correcta se conoce por
construcción y no depende de ninguna simulación. Ejecutar antes de usar el
módulo con datos reales:  python3 scripts/test_branch_counts.py
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from branch_counts import branch_counts, jackknife_fraction  # noqa: E402

fails = []


def check(name, got, want, tol=0):
    ok = (abs(got - want) <= tol) if isinstance(want, float) else (got == want)
    print(f"  {'ok  ' if ok else 'FALLA'} {name}: obtenido {got}, esperado {want}")
    if not ok:
        fails.append(name)


print("1. Caso mínimo: conjuntos de referencia bien separados")
#            sitio:  0  1  2  3  4  5  6  7
a = np.array([1, 1, 1, 0, 0, 0, 1, 0])
b = np.array([1, 1, 1, 1, 1, 1, 0, 0])
m = np.array([0, 0, 0, 0, 0, 0, 0, 0])
# AB = sitios 0,1,2 ; B = 3,4,5 ; A = 6 ; sitio 7 en ninguno
x = np.array([1, 1, 0, 1, 0, 0, 1, 1])
r = branch_counts(x, a, b, m)
check("n del conjunto AB", r["AB"]["n"], 3)
check("k del conjunto AB", r["AB"]["k"], 2)
check("f del conjunto AB", r["AB"]["f"], 2 / 3, 1e-12)
check("n del conjunto A", r["A"]["n"], 1)
check("f del conjunto A", r["A"]["f"], 1.0, 1e-12)
check("n del conjunto B", r["B"]["n"], 3)
check("f del conjunto B", r["B"]["f"], 1 / 3, 1e-12)

print("\n2. Un sitio donde M lleva el derivado queda excluido del conjunto de referencia")
a = np.array([1, 1]); b = np.array([1, 1]); m = np.array([0, 1]); x = np.array([1, 1])
r = branch_counts(x, a, b, m)
check("solo cuenta el sitio con M ancestral", r["AB"]["n"], 1)

print("\n3. Sin dato en X reduce el denominador, no lo cuenta como ancestral")
a = np.array([1, 1, 1, 1]); b = np.array([1, 1, 1, 1]); m = np.zeros(4, int)
x = np.array([1, 1, -1, -1])
r = branch_counts(x, a, b, m)
check("n con dos sitios sin dato", r["AB"]["n"], 2)
check("f no se diluye", r["AB"]["f"], 1.0, 1e-12)

print("\n4. Sin dato en A o B excluye el sitio por completo")
a = np.array([1, -1]); b = np.array([1, 1]); m = np.zeros(2, int); x = np.array([1, 1])
r = branch_counts(x, a, b, m)
check("n ignora el sitio sin dato en A", r["AB"]["n"], 1)

print("\n5. Proporción recuperada en una muestra grande con fracción fijada")
rng = np.random.default_rng(20261025)
N = 200_000
a = np.ones(N, int); b = np.ones(N, int); m = np.zeros(N, int)
for target in (0.91, 0.45, 0.20):
    x = (rng.random(N) < target).astype(int)
    r = branch_counts(x, a, b, m)
    check(f"f recupera {target}", round(r["AB"]["f"], 3), target, 0.004)

print("\n6. El error típico del jackknife concuerda con el binomial cuando los bloques son homogéneos")
N = 120_000
a = np.ones(N, int); b = np.ones(N, int); m = np.zeros(N, int)
x = (rng.random(N) < 0.45).astype(int)
blocks = np.repeat(np.arange(600), N // 600)
r = branch_counts(x, a, b, m, blocks=blocks)
f = r["AB"]["f"]
ee_binom = np.sqrt(f * (1 - f) / N)
ratio = r["AB"]["ee"] / ee_binom
check("cociente EE jackknife / EE binomial", round(ratio, 2), 1.0, 0.15)

print("\n7. El jackknife detecta estructura por bloques: EE mayor que el binomial")
#   mitad de los bloques con f alta, mitad con f baja -> el binomial subestima
blocks = np.repeat(np.arange(600), N // 600)
hi = blocks < 300
x = np.where(hi, (rng.random(N) < 0.70), (rng.random(N) < 0.20)).astype(int)
r = branch_counts(x, a, b, m, blocks=blocks)
f = r["AB"]["f"]
ee_binom = np.sqrt(f * (1 - f) / N)
mayor = r["AB"]["ee"] > 3 * ee_binom
print(f"  {'ok  ' if mayor else 'FALLA'} EE jackknife {r['AB']['ee']:.5f} frente a binomial {ee_binom:.5f}")
if not mayor:
    fails.append("estructura por bloques")

print("\n8. Entradas inválidas producen error")
for bad in ([2, 0], [0.5, 0]):
    try:
        branch_counts(np.array(bad), np.array([1, 1]), np.array([1, 1]), np.array([0, 0]))
        print(f"  FALLA no detecta la entrada {bad}")
        fails.append("validación")
    except ValueError:
        print(f"  ok   rechaza la entrada {bad}")
try:
    branch_counts(np.array([1]), np.array([1, 1]), np.array([1, 1]), np.array([0, 0]))
    print("  FALLA no detecta largos distintos")
    fails.append("validación de largo")
except ValueError:
    print("  ok   rechaza arrays de largo distinto")

print()
if fails:
    print(f"RESULTADO: {len(fails)} prueba(s) fallidas: {', '.join(fails)}")
    sys.exit(1)
print("RESULTADO: todas las pruebas pasan")
