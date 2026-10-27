#!/usr/bin/env python3
"""Pruebas de respuesta conocida para read_vcfs.py.

Construye cinco VCF diminutos que imitan el formato real de los ficheros de MPI-EVA, con
casos cuya respuesta se conoce por construcción, ejecuta el lector y comprueba la salida.

    python3 scripts/test_read_vcfs.py
"""
import gzip
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


ARCH_HDR = """##fileformat=VCFv4.1
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{name}
"""
MB_HDR = """##fileformat=VCFv4.1
##INFO=<ID=CAnc,Number=1,Type=String,Description="Ref-Chimp/Human ancestor base">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSS6004471
"""


def arch_rec(pos, ref, alt, gt, dp=50):
    return f"21\t{pos}\t.\t{ref}\t{alt}\t100\t.\t.\tGT:DP:A:C:G:T:PP:GQ\t{gt}:{dp}:0,0:0,0:0,0:0,0:0:99\n"


def mb_rec(pos, ref, alt, gt, info):
    return f"21\t{pos}\t.\t{ref}\t{alt}\t100\t.\t{info}\tGT:DP:GQ:PL\t{gt}:50:99:0,0,0\n"


def write(path, header, recs):
    with gzip.open(path, "wt") as fh:
        fh.write(header)
        fh.writelines(recs)


def build(d):
    """Casos. REF = A en todas las posiciones, CAnc = A salvo donde se indique.

    pos 100: todos 0/0            -> todos ancestral (0)
    pos 200: D5 y Vi 1/1 (alt G)  -> D5=1, Vi=1, D3=0, Chag8=0, Mbuti=0
    pos 300: sin CAnc             -> NO se emite
    pos 400: CAnc en minúscula    -> se emite, se trata igual
    pos 500: D3 ausente del fichero -> D3 = -1
    pos 600: Chag8 multialélico 1/2 -> Chag8 = -1
    pos 700: marcada RM           -> se emite salvo con --exclude-rm
    pos 800: D5 heterocigoto 0/1  -> en modo homozygous, D5 = -1
    pos 900: Mbuti derivado 1/1   -> Mbuti = 1
    """
    common = [(100, "0/0", "."), (200, "0/0", "."), (300, "0/0", "."), (400, "0/0", "."),
              (500, "0/0", "."), (600, "0/0", "."), (700, "0/0", "."), (800, "0/0", "."),
              (900, "0/0", ".")]

    def recs(overrides, skip=()):
        out = []
        for pos, gt, alt in common:
            if pos in skip:
                continue
            g, al = overrides.get(pos, (gt, alt))
            out.append(arch_rec(pos, "A", al, g))
        return out

    write(d / "D5_chr21.vcf.gz", ARCH_HDR.format(name="AltaiNeandertal"),
          recs({200: ("1/1", "G"), 800: ("0/1", "G")}))
    write(d / "Vi3319_chr21.vcf.gz", ARCH_HDR.format(name="Vindija33.19"),
          recs({200: ("1/1", "G")}))
    write(d / "D3_chr21.vcf.gz", ARCH_HDR.format(name="Denisova"),
          recs({}, skip=(500,)))
    write(d / "Chag8_chr21.vcf.gz", ARCH_HDR.format(name="Chagyrskaya-Phalanx"),
          recs({600: ("1/2", "G,T")}))

    mb = [
        mb_rec(100, "A", ".", "0/0", "AC=0;CAnc=A;Map20=1"),
        mb_rec(200, "A", ".", "0/0", "AC=0;CAnc=A;Map20=1"),
        mb_rec(300, "A", ".", "0/0", "AC=0;Map20=1"),                 # sin CAnc
        mb_rec(400, "A", ".", "0/0", "AC=0;CAnc=a;Map20=1"),          # minúscula
        mb_rec(500, "A", ".", "0/0", "AC=0;CAnc=A;Map20=1"),
        mb_rec(600, "A", ".", "0/0", "AC=0;CAnc=A;Map20=1"),
        mb_rec(700, "A", ".", "0/0", "AC=0;CAnc=A;RM;Map20=1"),       # repeat masked
        mb_rec(800, "A", ".", "0/0", "AC=0;CAnc=A;Map20=1"),
        mb_rec(900, "A", "G", "1/1", "AC=2;CAnc=A;Map20=1"),          # Mbuti derivado
    ]
    write(d / "Mbuti_chr21.vcf.gz", MB_HDR, mb)


def run(d, out, *extra):
    cmd = [sys.executable, str(HERE / "read_vcfs.py"), "--dir", str(d), "--out", str(out),
           "--mode", "homozygous", "--progress", "0", *extra]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout); print(r.stderr, file=sys.stderr)
        raise SystemExit("read_vcfs.py terminó con error")
    return np.load(out, allow_pickle=False)


with tempfile.TemporaryDirectory() as td:
    d = Path(td)
    build(d)

    print("A. Ejecución normal (modo homozygous)")
    z = run(d, d / "out.npz")
    pos = list(z["pos"])
    check("posiciones emitidas", pos, [100, 200, 400, 500, 600, 700, 800, 900])
    check("pos 300 sin CAnc queda fuera", 300 in pos, False)

    idx = {p: i for i, p in enumerate(pos)}

    def g(name, p):
        return int(z[name][idx[p]])

    print("\nB. Codificación ancestral / derivado")
    check("pos 100, D5 ancestral", g("D5", 100), 0)
    check("pos 200, D5 derivado", g("D5", 200), 1)
    check("pos 200, Vi derivado", g("Vi3319", 200), 1)
    check("pos 200, D3 ancestral", g("D3", 200), 0)
    check("pos 200, Chag8 ancestral", g("Chag8", 200), 0)
    check("pos 200, Mbuti ancestral", g("Mbuti", 200), 0)
    check("pos 900, Mbuti derivado", g("Mbuti", 900), 1)
    check("pos 400, CAnc minúscula se acepta", g("D5", 400), 0)

    print("\nC. Datos ausentes")
    check("pos 500, D3 ausente del fichero", g("D3", 500), -1)
    check("pos 500, D5 sí tiene dato", g("D5", 500), 0)
    check("pos 600, Chag8 multialélico", g("Chag8", 600), -1)
    check("pos 800, D5 heterocigoto en modo homozygous", g("D5", 800), -1)
    check("pos 800, Vi homocigoto sí cuenta", g("Vi3319", 800), 0)

    print("\nD. Filtro RM")
    z2 = run(d, d / "out_rm.npz", "--exclude-rm")
    check("pos 700 descartada con --exclude-rm", 700 in list(z2["pos"]), False)
    check("las demás siguen", list(z2["pos"]), [100, 200, 400, 500, 600, 800, 900])

    print("\nE. Filtro de profundidad")
    z3 = run(d, d / "out_dp.npz", "--min-dp", "100")
    check("con --min-dp 100 ningún arcaico tiene dato",
          int((z3["D5"] >= 0).sum()), 0)
    check("Mbuti no se ve afectado por --min-dp",
          int((z3["Mbuti"] >= 0).sum()), len(z3["pos"]))

    print("\nF. El modo pseudohaploid respeta la semilla")
    p1 = run(d, d / "p1.npz", "--mode", "pseudohaploid", "--seed", "7")
    p2 = run(d, d / "p2.npz", "--mode", "pseudohaploid", "--seed", "7")
    check("dos ejecuciones con la misma semilla coinciden",
          bool((p1["D5"] == p2["D5"]).all()), True)
    check("el heterocigoto de pos 800 ahora sí tiene dato",
          int(p1["D5"][list(p1["pos"]).index(800)]) in (0, 1), True)

    print("\nH. Modo carrier: el heterocigoto cuenta como derivado")
    zc = run(d, d / "carrier.npz", "--mode", "carrier")
    ic = {int(p): i for i, p in enumerate(zc["pos"])}
    check("pos 800, D5 heterocigoto (0/1) -> derivado",
          int(zc["D5"][ic[800]]), 1)
    check("pos 100, D5 homocigoto ancestral -> ancestral",
          int(zc["D5"][ic[100]]), 0)
    check("pos 200, D5 homocigoto derivado -> derivado",
          int(zc["D5"][ic[200]]), 1)
    check("pos 500, D3 ausente sigue siendo sin dato",
          int(zc["D3"][ic[500]]), -1)

    print("\nG. Encadenado con branch_counts")
    sys.path.insert(0, str(HERE))
    from branch_counts import branch_counts  # noqa: E402
    r = branch_counts(z["Chag8"], z["D5"], z["Vi3319"], z["Mbuti"])
    # conjunto AB = posiciones donde D5 y Vi son derivados y Mbuti ancestral -> solo pos 200
    check("conjunto AB tiene 1 sitio", r["AB"]["n"], 1)
    check("Chag8 es ancestral ahí, f_AB = 0", r["AB"]["f"], 0.0)

    print("\nI. Filtro --only-informative")
    zi = run(d, d / "inf.npz", "--only-informative")
    pi = [int(x) for x in zi["pos"]]
    # se emiten solo las posiciones con Mbuti ancestral y algun arcaico derivado: 200 y 800*
    # (800 es heterocigoto en D5 y en modo homozygous queda sin dato, asi que no cuenta)
    check("solo quedan posiciones informativas", pi, [200])
    check("la base ancestral se guarda", chr(int(zi["canc"][0])), "A")
    check("la base derivada se guarda", chr(int(zi["der"][0])), "G")
    check("el cromosoma se guarda", str(zi["chrom"][0]), "21")

print()
if fails:
    print(f"RESULTADO: {len(fails)} prueba(s) fallidas: {', '.join(fails)}")
    sys.exit(1)
print("RESULTADO: todas las pruebas pasan")
