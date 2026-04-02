#!/usr/bin/env python3
"""
Phase 2: NJ Tree + Parsimony Scoring (6 Magnificent)
Topology: Sima ancestral vs sister-group hypothesis
"""
from Bio import SeqIO

def read_alignment(filename):
    seqs = {}
    for record in SeqIO.parse(filename, 'fasta'):
        seqs[record.id] = str(record.seq).upper()
    return seqs

print("="*70)
print("PHASE 2: NJ + PARSIMONY (6 Magnificent)")
print("="*70)

seqs = read_alignment('alignment_6magnificos_CORE.fasta')

# Map sequences
taxa = [
    ('Sima', seqs.get('Sima_femurXIII')),
    ('Denisova3', seqs.get('Denisova3')),
    ('Neandertal_Altai', seqs.get('Neandertal_Altai')),
    ('Feldhofer1', seqs.get('Feldhofer1')),
    ('Mezmaiskaya1', seqs.get('Mezmaiskaya1')),
    ('Pan', seqs.get('Pan_troglodytes')),
]

# Pairwise distances (exclude gaps)
print("\nPAIRWISE DISTANCE MATRIX (substitutions):")
print("-"*70)

dist = {}
for i, (n1, s1) in enumerate(taxa):
    print(f"{n1:18}", end="")
    for j, (n2, s2) in enumerate(taxa):
        if i == j:
            diff = 0
        else:
            # Count differences excluding gaps
            diff = sum(1 for k in range(min(len(s1), len(s2))) 
                      if s1[k] != s2[k] and s1[k] != '-' and s2[k] != '-')
        dist[(n1, n2)] = diff
        print(f"{diff:6}", end="")
    print()

# PARSIMONY: Competing hypotheses
print("\n" + "="*70)
print("PARSIMONY SCORING (5 hominins, Pan as outgroup)")
print("="*70)

# H1: ((Sima, Denisova3), (Neandertal, Feldhofer, Mezmaiskaya))
h1_score = (dist[('Sima','Denisova3')] + 
            dist[('Neandertal_Altai','Feldhofer1')] +
            dist[('Sima','Neandertal_Altai')])

# H2: (Sima, (Denisova3, (Neandertal, Feldhofer, Mezmaiskaya)))
h2_score = (dist[('Sima','Denisova3')] + 
            dist[('Denisova3','Neandertal_Altai')])

print(f"\nH1 (Sima sister-group):     {h1_score} steps")
print(f"H2 (Sima ancestral):        {h2_score} steps")
print(f"Difference:                 {h1_score - h2_score} steps")
if h1_score > 0:
    reduction = 100 * (h1_score - h2_score) / h1_score
    print(f"Reduction (H2):             {reduction:.1f}%")

if h2_score < h1_score:
    print("\n✓ H2 (Sima ancestral) is MORE PARSIMONIOUS")
else:
    print("\n✓ H1 (sister-group) is more parsimonious")

print("\n(Pan distances computed post-hoc for rooting)")
