#!/usr/bin/env python3
"""Comprueba que cada npz del recorrido genómico cubre su cromosoma entero.

Un VCF truncado produce un npz que parece válido pero solo abarca el principio del
cromosoma. Esto lo detecta comparando la posición máxima con la longitud conocida.

    python3 scripts/check_genome.py results/genome
"""
import sys
from pathlib import Path

import numpy as np

# longitudes hg19 (GRCh37)
LEN = {1: 249250621, 2: 243199373, 3: 198022430, 4: 191154276, 5: 180915260,
       6: 171115067, 7: 159138663, 8: 146364022, 9: 141213431, 10: 135534747,
       11: 135006516, 12: 133851895, 13: 115169878, 14: 107349540, 15: 102531392,
       16: 90354753, 17: 81195210, 18: 78077248, 19: 59128983, 20: 63025520,
       21: 48129895, 22: 51304566}

d = Path(sys.argv[1] if len(sys.argv) > 1 else "results/genome")
print(f"{'cr':>3s} {'sitios':>11s} {'pos max':>12s} {'long cr':>12s} {'cubre':>7s}  estado")
malos, total = [], 0
for c in sorted(LEN):
    f = d / f"chr{c}_inf.npz"
    if not f.exists():
        continue
    z = np.load(f, allow_pickle=False)
    n = len(z["pos"])
    total += n
    mx = int(z["pos"].max()) if n else 0
    frac = mx / LEN[c]
    # los VCF empiezan tras el centrómero en los acrocéntricos, pero el final debe
    # quedar cerca del extremo del cromosoma en todos los casos
    bad = frac < 0.90
    if bad:
        malos.append(c)
    print(f"{c:3d} {n:11,} {mx:12,} {LEN[c]:12,} {frac:6.1%}  "
          f"{'TRUNCADO' if bad else 'ok'}")

print(f"\nTotal: {total:,} sitios informativos")
if malos:
    print(f"\nTRUNCADOS: {malos}")
    print("Bórralos y relanza el recorrido, que los rehará:")
    print("  rm " + " ".join(f"results/genome/chr{c}_inf.npz" for c in malos))
    sys.exit(1)
print("Todos los cromosomas presentes cubren su longitud completa.")
