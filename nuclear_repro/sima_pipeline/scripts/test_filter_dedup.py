#!/usr/bin/env python3
"""Pruebas de respuesta conocida para filter_dedup.py.

Cada caso se construye a mano con un SAM mínimo, de modo que la respuesta correcta se
conoce por construcción.

    python3 scripts/test_filter_dedup.py
"""
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).parent
fails = []


def check(name, got, want):
    ok = got == want
    print(f"  {'ok  ' if ok else 'FALLA'} {name}: obtenido {got}, esperado {want}")
    if not ok:
        fails.append(name)


HDR = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:21\tLN:48129895\n"


def rec(nombre, flag, pos, mapq, cigar, seq):
    q = "I" * len(seq)
    return f"{nombre}\t{flag}\t21\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{seq}\t{q}\n"


def run(sam, *args):
    r = subprocess.run([sys.executable, str(HERE / "filter_dedup.py"), *args],
                       input=sam, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stderr, file=sys.stderr)
        raise SystemExit("filter_dedup.py terminó con error")
    nombres = [l.split("\t")[0] for l in r.stdout.splitlines() if not l.startswith("@")]
    cab = [l for l in r.stdout.splitlines() if l.startswith("@")]
    return nombres, cab, r.stderr


S40 = "A" * 40
S34 = "A" * 34

print("1. La cabecera se conserva intacta")
nombres, cab, _ = run(HDR + rec("r1", 0, 1000, 60, "40M", S40))
check("líneas de cabecera", len(cab), 2)
check("lecturas conservadas", nombres, ["r1"])

print("\n2. Filtros básicos")
sam = (HDR
       + rec("buena", 0, 1000, 60, "40M", S40)
       + rec("sin_mapear", 4, 0, 0, "*", S40)
       + rec("mapq_bajo", 0, 2000, 24, "40M", S40)
       + rec("mapq_justo", 0, 2100, 25, "40M", S40)
       + rec("corta", 0, 3000, 60, "34M", S34)
       + rec("secundaria", 256, 4000, 60, "40M", S40)
       + rec("suplementaria", 2048, 5000, 60, "40M", S40))
nombres, _, err = run(sam)
check("pasan las esperadas", nombres, ["buena", "mapq_justo"])
check("mapq 25 pasa (el umbral es >=)", "mapq_justo" in nombres, True)
check("mapq 24 no pasa", "mapq_bajo" in nombres, False)
check("longitud 34 no pasa", "corta" in nombres, False)

print("\n3. Duplicados: misma posición, mismo fin, misma hebra")
sam = (HDR
       + rec("original", 0, 1000, 60, "40M", S40)
       + rec("copia", 0, 1000, 60, "40M", S40)
       + rec("otra_copia", 0, 1000, 60, "40M", S40))
nombres, _, _ = run(sam)
check("solo queda una", nombres, ["original"])

print("\n4. NO son duplicados si difiere algo de la clave")
sam = (HDR
       + rec("ref", 0, 1000, 60, "40M", S40)
       + rec("otro_inicio", 0, 1001, 60, "40M", S40)
       + rec("otra_hebra", 16, 1000, 60, "40M", S40)
       + rec("otro_fin", 0, 1000, 60, "45M", "A" * 45))
nombres, _, _ = run(sam)
check("se conservan las cuatro", sorted(nombres),
      ["otra_hebra", "otro_fin", "otro_inicio", "ref"])

print("\n5. El fin se calcula con el CIGAR, no con la longitud de la secuencia")
# 40M y 20M5D20M ocupan 40 y 45 pb sobre la referencia aunque ambas tengan 40 pb de
# secuencia: no deben considerarse duplicadas.
sam = (HDR
       + rec("sin_delecion", 0, 1000, 60, "40M", S40)
       + rec("con_delecion", 0, 1000, 60, "20M5D20M", S40))
nombres, _, _ = run(sam)
check("la deleción cambia el fin", sorted(nombres), ["con_delecion", "sin_delecion"])

# Insercion: 20M5I15M consume 35 pb de referencia con 40 de secuencia
sam = (HDR
       + rec("sin_insercion", 0, 2000, 60, "35M", "A" * 35)
       + rec("con_insercion", 0, 2000, 60, "20M5I15M", S40))
nombres, _, _ = run(sam)
check("la inserción da el mismo fin: son duplicadas", nombres, ["sin_insercion"])

print("\n6. El recorte blando (S) no consume referencia")
sam = (HDR
       + rec("recortada", 0, 3000, 60, "5S40M", "A" * 45)
       + rec("entera", 0, 3000, 60, "40M", S40))
nombres, _, _ = run(sam)
check("ambas acaban en el mismo sitio: duplicadas", nombres, ["recortada"])

print("\n7. --no-dedup conserva los duplicados")
sam = HDR + rec("a", 0, 1000, 60, "40M", S40) + rec("b", 0, 1000, 60, "40M", S40)
nombres, _, _ = run(sam, "--no-dedup")
check("se conservan las dos", nombres, ["a", "b"])

print("\n8. Umbrales configurables")
sam = HDR + rec("m30", 0, 1000, 60, "40M", S40) + rec("m30b", 0, 2000, 30, "40M", S40)
nombres, _, _ = run(sam, "--min-mapq", "37")
check("con --min-mapq 37 solo pasa la de 60", nombres, ["m30"])
sam = HDR + rec("l40", 0, 1000, 60, "40M", S40)
nombres, _, _ = run(sam, "--min-len", "45")
check("con --min-len 45 no pasa ninguna", nombres, [])

print("\n9. Las estadísticas cuadran")
sam = (HDR
       + rec("a", 0, 1000, 60, "40M", S40)
       + rec("b", 0, 1000, 60, "40M", S40)
       + rec("c", 4, 0, 0, "*", S40)
       + rec("d", 0, 2000, 10, "40M", S40))
nombres, _, err = run(sam)
def val(etq):
    for l in err.splitlines():
        if etq in l:
            return int(l.split(":")[1].strip().replace(",", ""))
    return None
check("leídas", val("lecturas leídas"), 4)
check("sin mapear", val("sin mapear"), 1)
check("calidad insuficiente", val("calidad insuficiente"), 1)
check("duplicadas", val("duplicadas"), 1)
check("conservadas", val("CONSERVADAS"), 1)

print()
if fails:
    print(f"RESULTADO: {len(fails)} prueba(s) fallidas: {', '.join(fails)}")
    sys.exit(1)
print("RESULTADO: todas las pruebas pasan")
