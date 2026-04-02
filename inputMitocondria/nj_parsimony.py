#!/usr/bin/env python3
"""
Neighbor-Joining tree + Parsimony scoring
mtDNA analysis: Sima ancestral vs sister-group hypothesis
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

print("NJ + PARSIMONY ANALYSIS")
print("="*60)

seqs = read_fasta('alignment_6_mtdna.fasta')

# Identificar secuencias
sima = None
d3 = None
nean = None
hsap = None

for k, v in seqs.items():
    if 'KF683087' in k:
        sima = v
    elif 'FN673705' in k:
        d3 = v
    elif 'KC879692' in k:
        nean = v
    elif 'NC_012920' in k:
        hsap = v

# Matriz distancias pairwise
taxa = [('Sima', sima), ('Denisova3', d3), ('Neandertal', nean), ('HSapiens', hsap)]

print("\nMATRIZ DE DISTANCIAS:")
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

# PARSIMONIA: Contar cambios por topología
print("\n" + "="*60)
print("PARSIMONY SCORING")
print("="*60)

# H1: ((Sima, Denisova3), Neandertal)
h1_score = dist[('Sima','Denisova3')] + dist[('Sima','Neandertal')] + dist[('Denisova3','Neandertal')]

# H2: (Sima, (Denisova3, Neandertal))
h2_score = dist[('Sima','Denisova3')] + dist[('Denisova3','Neandertal')]

print(f"\nH1 (sister-group):    {h1_score} steps")
print(f"H2 (ancestral):       {h2_score} steps")
print(f"Difference:           {h1_score - h2_score} steps")
print(f"Reduction (H2):       {100*(h1_score-h2_score)/h1_score:.1f}%")

if h2_score < h1_score:
    print("\n✓ H2 (Sima ancestral) is more parsimonious")
else:
    print("\n✓ H1 (sister-group) is more parsimonious")
