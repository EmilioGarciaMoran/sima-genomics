#!/usr/bin/env python3
"""Pruebas de respuesta conocida para damage_filter.py y fasta.py.

Construye un FASTA diminuto con su .fai hecho a mano y lecturas sintéticas con el daño
puesto a propósito, de modo que la respuesta correcta se conoce por construcción.

    python3 scripts/test_damage_filter.py
"""
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
fails = []


def check(name, got, want):
    ok = got == want
    print(f"  {'ok  ' if ok else 'FALLA'} {name}: obtenido {got}, esperado {want}")
    if not ok:
        fails.append(name)


# --- FASTA de prueba: 120 pb con un patrón conocido ----------------------
REF = ("CATGCATGCATGCATGCATG" * 6)          # 120 pb, periodo 4: C A T G
assert len(REF) == 120


def escribir_fasta(d, bpl=60):
    fa = d / "t.fasta"
    with open(fa, "w") as f:
        f.write(">21 prueba\n")
        for i in range(0, len(REF), bpl):
            f.write(REF[i:i+bpl] + "\n")
    # .fai: nombre, longitud, desplazamiento, bases por línea, bytes por línea
    off = len(">21 prueba\n")
    with open(str(fa) + ".fai", "w") as f:
        f.write(f"21\t{len(REF)}\t{off}\t{bpl}\t{bpl+1}\n")
    return fa


HDR = "@HD\tVN:1.6\n@SQ\tSN:21\tLN:120\n"


def rec(nombre, flag, pos, cigar, seq):
    return f"{nombre}\t{flag}\t21\t{pos}\t60\t{cigar}\t*\t0\t0\t{seq}\t{'I'*len(seq)}\n"


def run(fa, sam, *args):
    r = subprocess.run([sys.executable, str(HERE / "damage_filter.py"), "--ref", str(fa), *args],
                       input=sam, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stderr, file=sys.stderr)
        raise SystemExit("damage_filter.py terminó con error")
    datos = [l.split("\t") for l in r.stdout.splitlines() if not l.startswith("@")]
    return {d[0]: d[9] for d in datos}, r.stderr


with tempfile.TemporaryDirectory() as td:
    d = Path(td)
    fa = escribir_fasta(d)

    print("1. El lector de FASTA devuelve las bases correctas")
    from fasta import Fasta
    f = Fasta(fa)
    check("primeras 8 bases", f.fetch("21", 0, 8), REF[:8])
    check("a caballo entre dos líneas", f.fetch("21", 55, 65), REF[55:65])
    check("últimas 5", f.fetch("21", 115, 120), REF[115:])
    check("más allá del final se recorta", f.fetch("21", 115, 200), REF[115:])
    check("cromosoma inexistente", f.fetch("X", 0, 10), "")
    check("longitud declarada", f.length("21"), 120)
    f.close()

    # El patrón CATG repetido: posicion 0-based 0=C, 1=A, 2=T, 3=G, 4=C, ...
    # Una lectura de 40 pb desde la posicion 0 (1-based: 1) empieza en C y acaba en...
    # indice 39 -> 39 % 4 = 3 -> G.  Perfecto: C al principio, G al final.
    SANA = REF[0:40]

    print("\n2. Cadena sencilla, hebra directa: se exige C->T en algún extremo")
    dañada_izq = "T" + SANA[1:]                 # C->T en la primera posicion
    dañada_der = SANA[:39] + "A"                # G->A al final: NO es daño en hebra directa
    sam = (HDR
           + rec("sana", 0, 1, "40M", SANA)
           + rec("dano_izq", 0, 1, "40M", dañada_izq)
           + rec("ga_derecha", 0, 1, "40M", dañada_der))
    seqs, err = run(fa, sam, "--mask", "0")
    check("la sana se descarta", "sana" in seqs, False)
    check("la dañada se conserva", "dano_izq" in seqs, True)
    check("G->A no cuenta en hebra directa", "ga_derecha" in seqs, False)

    print("\n3. Cadena sencilla, hebra inversa: se exige G->A")
    sam = (HDR
           + rec("inv_sana", 16, 1, "40M", SANA)
           + rec("inv_ga", 16, 1, "40M", dañada_der)
           + rec("inv_ct", 16, 1, "40M", dañada_izq))
    seqs, _ = run(fa, sam, "--mask", "0")
    check("G->A se conserva en hebra inversa", "inv_ga" in seqs, True)
    check("C->T no cuenta en hebra inversa", "inv_ct" in seqs, False)
    check("la sana se descarta", "inv_sana" in seqs, False)

    print("\n4. Doble cadena: C->T a la izquierda, G->A a la derecha, sin mirar la hebra")
    sam = (HDR
           + rec("ds_ct", 0, 1, "40M", dañada_izq)
           + rec("ds_ga", 0, 1, "40M", dañada_der)
           + rec("ds_inv_ct", 16, 1, "40M", dañada_izq))
    seqs, _ = run(fa, sam, "--library", "ds", "--mask", "0")
    check("C->T izquierda vale", "ds_ct" in seqs, True)
    check("G->A derecha vale", "ds_ga" in seqs, True)
    check("en hebra inversa tambien vale C->T", "ds_inv_ct" in seqs, True)

    print("\n5. Enmascarado: T a N en los 3 primeros lugares (hebra directa)")
    # SANA = C A T G C A T G ...  la T esta en el indice 2, dentro de los 3 primeros
    sam = HDR + rec("r", 0, 1, "40M", dañada_izq)
    seqs, _ = run(fa, sam, "--mask", "3")
    s = seqs["r"]
    check("posicion 0 (T de daño) enmascarada", s[0], "N")
    check("posicion 1 (A) intacta", s[1], "A")
    check("posicion 2 (T) enmascarada", s[2], "N")
    check("posicion 3 (G) intacta", s[3], "G")
    check("la longitud no cambia", len(s), 40)

    print("\n6. Enmascarado en el extremo derecho (hebra directa: tambien T)")
    # indices 37,38,39 -> 37%4=1 A, 38%4=2 T, 39%4=3 G
    check("posicion 38 (T) enmascarada", s[38], "N")
    check("posicion 37 (A) intacta", s[37], "A")
    check("posicion 39 (G) intacta", s[39], "G")

    print("\n7. Hebra inversa: se enmascaran las A, no las T")
    sam = HDR + rec("r", 16, 1, "40M", dañada_der)
    seqs, _ = run(fa, sam, "--mask", "3")
    s = seqs["r"]
    check("posicion 1 (A) enmascarada", s[1], "N")
    check("posicion 2 (T) intacta", s[2], "T")
    check("posicion 39 (A de daño) enmascarada", s[39], "N")

    print("\n8. El recorte blando no cuenta como posición alineada")
    # 5S35M: las 5 primeras bases de SEQ no se alinean. La alineacion empieza en la
    # posicion 1 de la referencia, asi que la base alineada 0 corresponde a REF[0]=C.
    seq = "GGGGG" + "T" + SANA[1:35]
    sam = HDR + rec("recortada", 0, 1, "5S35M", seq)
    seqs, _ = run(fa, sam, "--mask", "3")
    check("se conserva (el daño está en la primera base alineada)", "recortada" in seqs, True)
    s = seqs["recortada"]
    check("el recorte no se enmascara", s[:5], "GGGGG")
    check("la primera alineada sí", s[5], "N")

    print("\n8 bis. Ventana de deteccion del daño")
    # SANA = C A T G C A T G...  el indice 4 es una C: ponerla a T es daño en la posicion 5
    dano_pos5 = SANA[:4] + "T" + SANA[5:]
    sam = HDR + rec("w", 0, 1, "40M", dano_pos5)
    seqs, _ = run(fa, sam, "--mask", "0", "--damage-window", "1")
    check("con ventana 1 no se detecta", "w" in seqs, False)
    seqs, _ = run(fa, sam, "--mask", "0", "--damage-window", "5")
    check("con ventana 5 si se detecta", "w" in seqs, True)
    # y en el extremo derecho: indice 36 es C (36%4=0)
    dano_der5 = SANA[:36] + "T" + SANA[37:]
    sam = HDR + rec("w2", 0, 1, "40M", dano_der5)
    seqs, _ = run(fa, sam, "--mask", "0", "--damage-window", "1")
    check("extremo derecho, ventana 1: no", "w2" in seqs, False)
    seqs, _ = run(fa, sam, "--mask", "0", "--damage-window", "4")
    check("extremo derecho, ventana 4: si", "w2" in seqs, True)

    print("\n9. --no-filter conserva todo y sigue enmascarando")
    sam = HDR + rec("sana", 0, 1, "40M", SANA) + rec("dano", 0, 1, "40M", dañada_izq)
    seqs, _ = run(fa, sam, "--no-filter", "--mask", "3")
    check("se conservan las dos", sorted(seqs), ["dano", "sana"])
    check("la sana también se enmascara", seqs["sana"][2], "N")

    print("\n10. Las estadísticas cuadran")
    sam = (HDR
           + rec("a", 0, 1, "40M", dañada_izq)
           + rec("b", 0, 1, "40M", SANA)
           + rec("c", 0, 1, "40M", dañada_izq))
    _, err = run(fa, sam, "--mask", "3")
    def val(etq):
        for l in err.splitlines():
            if etq in l:
                return int(l.split(":")[1].strip().replace(",", "").replace("%", ""))
        return None
    check("leídas", val("lecturas leídas"), 3)
    check("sin daño", val("sin daño terminal"), 1)
    check("conservadas", val("CONSERVADAS"), 2)
    check("bases alineadas", val("bases alineadas"), 80)

print()
if fails:
    print(f"RESULTADO: {len(fails)} prueba(s) fallidas: {', '.join(fails)}")
    sys.exit(1)
print("RESULTADO: todas las pruebas pasan")
