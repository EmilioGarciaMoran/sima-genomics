#!/usr/bin/env bash
# Paso 2: alineamiento múltiple con MAFFT y versión sin D-loop.
#
#   bash scripts/02_align.sh
#
# Entrada : sequences/all_raw.fasta
# Salidas : alignments/aln_full.fasta     genoma completo (MAFFT G-INS-i)
#           alignments/aln_coding.fasta   sin D-loop (posiciones 600-16000 de Vindija)
#           alignments/ALIGNMENT_REPORT.txt  versión, comando, estadísticas y comparación con el alineamiento legacy
# Ejecutar desde mtDNA_repro/ con el entorno mtphylo activo.
set -euo pipefail
IN=sequences/all_raw.fasta
OUT=alignments
LEGACY=../inputMitocondria/phase3a_results/alignment_6magnificos.phy
mkdir -p "$OUT"

MAFFT_CMD="mafft --globalpair --maxiterate 1000 --thread 1 $IN"
$MAFFT_CMD 2> /dev/null \
  | awk '/^>/{print $1; next}{print toupper($0)}' > "$OUT/aln_full.fasta"

python3 -I - "$OUT" "$LEGACY" "$MAFFT_CMD" "$(mafft --version 2>&1 | head -1)" <<'PY'
import sys
out, legacy, cmd, ver = sys.argv[1:]

def read_fa(p):
    s = {}
    for b in open(p).read().split('>')[1:]:
        h, q = b.split('\n', 1)
        s[h.split()[0]] = q.replace('\n', '').upper()
    return s

def stats(s):
    L = len(next(iter(s.values())))
    cols = list(zip(*s.values()))
    gap_cols = sum('-' in c for c in cols)
    var = sum(len({x for x in c if x in 'ACGT'}) > 1 for c in cols)
    pis = 0
    for c in cols:
        cnt = {}
        for x in c:
            if x in 'ACGT': cnt[x] = cnt.get(x, 0) + 1
        if sum(v >= 2 for v in cnt.values()) >= 2: pis += 1
    return L, gap_cols, var, pis

full = read_fa(f"{out}/aln_full.fasta")

# Sin D-loop: el D-loop en coordenadas rCRS ocupa 16024-16569 y 1-576.
# Se conservan las columnas en las que Vindija (numeración casi idéntica a rCRS) está entre 600 y 16000.
ref, pos, keep = full['Neandertal'], 0, []
for i, c in enumerate(ref):
    if c != '-': pos += 1
    if 600 <= pos <= 16000: keep.append(i)
coding = {n: ''.join(q[i] for i in keep) for n, q in full.items()}
with open(f"{out}/aln_coding.fasta", "w") as f:
    for n, q in coding.items(): f.write(f">{n}\n{q}\n")

# Comparación con el alineamiento legacy: ¿mismos pares homólogos?
leg = {}
for line in open(legacy).read().splitlines()[1:]:
    if line.strip(): leg[line[:10].strip()] = line[10:].upper()

def residue_columns(aln):
    # para cada columna, la tupla (índice del residuo en cada secuencia o None si es hueco)
    idx = {n: -1 for n in aln}
    cols = []
    for j in range(len(next(iter(aln.values())))):
        t = []
        for n in sorted(aln):
            if aln[n][j] != '-':
                idx[n] += 1; t.append(idx[n])
            else: t.append(None)
        cols.append(tuple(t))
    return set(cols)

a, b = residue_columns(full), residue_columns(leg)
shared = len(a & b)

with open(f"{out}/ALIGNMENT_REPORT.txt", "w") as f:
    f.write(f"MAFFT: {ver}\nComando: {cmd}\n\n")
    for name, s in (("aln_full", full), ("aln_coding", coding), ("legacy phase3a", leg)):
        L, g, v, p = stats(s)
        f.write(f"{name:15s} columnas={L:6d} columnas_con_hueco={g:5d} variables={v:5d} parsimonia_informativas={p:4d}\n")
    f.write(f"\nColumnas idénticas (mismos residuos homólogos) entre aln_full y legacy: {shared} de {len(a)} (aln_full) / {len(b)} (legacy)\n")
print(open(f"{out}/ALIGNMENT_REPORT.txt").read())
PY
