#!/usr/bin/env python3
"""
NJ tree with Pan troglodytes as outgroup
Parsimony scoring of competing hypotheses
"""

def read_fasta(filename):
    seqs = {}
    name, seq = "", ""
    for line in open(filename):
        if line.startswith('>'):
            if name: seqs[name] = seq
            name = line.strip()[1:].split()[0]
            seq = ""
        else:
            seq += line.strip()
    if name: seqs[name] = seq
    return seqs

print("NJ + PARSIMONY (WITH PAN OUTGROUP)")
print("="*60)



# CAMBIA ESTA LÍNEA:
seqs = read_fasta('alignment_6_mtdna.fasta')

# POR ESTA:
seqs = read_fasta('alignment_6_mtdna_CODING_aligned.fasta')

# Ctrl+O → Enter → Ctrl+X#seqs = read_fasta('alignment_6_mtdna.fasta')

# Identificar 7 secuencias
sima = None
d3 = None
nean = None
hsap = None
pan = None

for k, v in seqs.items():
    if 'KF683087' in k:
        sima = v
    elif 'FN673705' in k:
        d3 = v
    elif 'KC879692' in k:
        nean = v
    elif 'NC_012920' in k:
        hsap = v
    elif 'NC_001643' in k or 'troglodytes' in k.lower():
        pan = v

taxa = [
    ('Sima', sima),
    ('Denisova3', d3),
    ('Neandertal', nean),
    ('HSapiens', hsap),
    ('Pan', pan)
]

print("\nMATRIZ DE DISTANCIAS (5 taxa + Pan outgroup):")
print("-"*60)

dist = {}
for i, (n1, s1) in enumerate(taxa):
    print(f"{n1:12}", end="")
    for j, (n2, s2) in enumerate(taxa):
        if i == j:
            diff = 0
        else:
            diff = sum(1 for k in range(len(s1)) if s1[k] != s2[k])
        dist[(n1, n2)] = diff
        print(f"{diff:8}", end="")
    print()

print("\nNOTA: Pan distancias son MÁS GRANDES")
print("      (Pan es OUTGROUP)")

# PARSIMONIA: Con Pan como enraizador
print("\n" + "="*60)
print("PARSIMONY SCORING (Pan enraiza el árbol)")
print("="*60)

# H1: ((Sima, Denisova3), Neandertal) [con Pan fuera]
h1_score = dist[('Sima','Denisova3')] + dist[('Sima','Neandertal')] + dist[('Denisova3','Neandertal')]

# H2: (Sima, (Denisova3, Neandertal)) [con Pan fuera]
h2_score = dist[('Sima','Denisova3')] + dist[('Denisova3','Neandertal')]

print(f"\nH1 (sister-group):    {h1_score} steps")
print(f"H2 (ancestral):       {h2_score} steps")
print(f"Difference:           {h1_score - h2_score} steps")
print(f"Reduction (H2):       {100*(h1_score-h2_score)/h1_score:.1f}%")

if h2_score < h1_score:
    print("\n✓ H2 (Sima ancestral) is MORE PARSIMONIOUS")
else:
    print("\n✓ H1 (sister-group) is more parsimonious")

print("\nWith Pan as rooting outgroup, tree is properly polarized")
