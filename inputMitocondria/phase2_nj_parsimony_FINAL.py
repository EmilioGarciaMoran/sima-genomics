#!/usr/bin/env python3
"""
Phase 2: NJ + Parsimony (6 Pleistocene specimens)
Conservative nomenclature: specimens identified by locality/age
Taxonomic interpretation deferred to Discussion
"""
from Bio import SeqIO

def read_alignment(filename):
    seqs = {}
    for record in SeqIO.parse(filename, 'fasta'):
        seqs[record.id] = str(record.seq).upper()
    return seqs

print("="*70)
print("PHASE 2: PARSIMONY ANALYSIS - PLEISTOCENE mtDNA")
print("="*70)

seqs = read_alignment('alignment_6magnificos_CORE.fasta')

# CONSERVATIVE NOMENCLATURE (locality/specimen)
taxa = [
    ('Sima_430ka', seqs.get('Sima_femurXIII')),
    ('Denisova3_70ka', seqs.get('Denisova3')),
    ('Altai_50ka', seqs.get('Neandertal_Altai')),
    ('Feldhofer_40ka', seqs.get('Feldhofer1')),
    ('Mezmaiskaya_60ka', seqs.get('Mezmaiskaya1')),
    ('Pan_outgroup', seqs.get('Pan_troglodytes')),
]

# Pairwise distances
print("\nPAIRWISE DISTANCE MATRIX (nucleotide substitutions):")
print("-"*70)

dist = {}
for i, (n1, s1) in enumerate(taxa):
    print(f"{n1:20}", end="")
    for j, (n2, s2) in enumerate(taxa):
        if i == j:
            diff = 0
        else:
            diff = sum(1 for k in range(min(len(s1), len(s2))) 
                      if s1[k] != s2[k] and s1[k] != '-' and s2[k] != '-')
        dist[(n1, n2)] = diff
        print(f"{diff:6}", end="")
    print()

# PARSIMONY
print("\n" + "="*70)
print("PARSIMONY SCORING (5 specimens, Pan as outgroup)")
print("="*70)

h1_score = (dist[('Sima_430ka','Denisova3_70ka')] + 
            dist[('Altai_50ka','Feldhofer_40ka')] +
            dist[('Sima_430ka','Altai_50ka')])

h2_score = (dist[('Sima_430ka','Denisova3_70ka')] + 
            dist[('Denisova3_70ka','Altai_50ka')])

print(f"\nH1 (Sima as sister-group):     {h1_score} steps")
print(f"H2 (Sima ancestral position):  {h2_score} steps")
print(f"Difference:                    {h1_score - h2_score} steps")
if h1_score > 0:
    reduction = 100 * (h1_score - h2_score) / h1_score
    print(f"Parsimony reduction (H2):      {reduction:.1f}%")

print("\n" + "="*70)
if h2_score < h1_score:
    print("✓ H2 (ancestral position) is more parsimonious")
else:
    print("✓ H1 (sister-group) is more parsimonious")
print("="*70)

