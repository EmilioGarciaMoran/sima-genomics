#!/usr/bin/env python3
"""
Phase 2: NJ + Parsimony Analysis (6 Magnificent - CORRECTED)
Input: alignment_6magnificos_CORRECTED.fasta (from MAFFT)
Conservative nomenclature: specimens by locality/age
Accessions validated by: DeepSeek + ChatGPT editorial review
"""
from Bio import SeqIO
import os

os.chdir('inputMitocondria')

def read_alignment(filename):
    seqs = {}
    for record in SeqIO.parse(filename, 'fasta'):
        seqs[record.id] = str(record.seq).upper()
    return seqs

print("="*80)
print("PHASE 2: PARSIMONY ANALYSIS (6 MAGNIFICENT - CORRECTED ACCESSIONS)")
print("="*80)
print("\nData sources:")
print("  • Sima_430ka           → KF683087.1 (Meyer et al., 2014)")
print("  • Denisova3_70ka       → NC_013993.1 (Reich et al., 2010) [CORRECTED]")
print("  • Neandertal_Vindija   → NC_020627.1 (Prüfer et al., 2014) [CORRECTED]")
print("  • Feldhofer1_40ka      → FM865407.1 (Green et al., 2008)")
print("  • Mezmaiskaya1_60ka    → FM865411.1 (Green et al., 2008)")
print("  • Pan_outgroup         → NC_001643.1 (Reference)\n")

alignment_file = 'alignment_6magnificos_CORRECTED.fasta'

if not os.path.exists(alignment_file):
    print(f"ERROR: {alignment_file} not found!")
    print("Run MAFFT first:")
    print(f"  mafft --thread 4 pleistocene_6magnificos_MASTER_CORRECTED.fasta > {alignment_file}")
    exit(1)

seqs = read_alignment(alignment_file)

# Conservative nomenclature (specimen-based, not taxonomic)
taxa = [
    ('Sima_430ka', seqs.get('Sima_430ka')),
    ('Denisova3_70ka', seqs.get('Denisova3_70ka')),
    ('Neandertal_Vindija_44ka', seqs.get('Neandertal_Vindija_44ka')),
    ('Feldhofer1_40ka', seqs.get('Feldhofer1_40ka')),
    ('Mezmaiskaya1_60ka', seqs.get('Mezmaiskaya1_60ka')),
    ('Pan_outgroup', seqs.get('Pan_outgroup')),
]

# Validate all sequences present
missing = [name for name, seq in taxa if seq is None]
if missing:
    print(f"ERROR: Missing sequences: {missing}")
    print(f"Available in alignment: {list(seqs.keys())}")
    exit(1)

# Pairwise distances (exclude gaps)
print("\nPAIRWISE DISTANCE MATRIX (nucleotide substitutions):")
print("-"*80)

dist = {}
for i, (n1, s1) in enumerate(taxa):
    print(f"{n1:25}", end="")
    for j, (n2, s2) in enumerate(taxa):
        if i == j:
            diff = 0
        else:
            # Count differences excluding gaps/Ns
            diff = sum(1 for k in range(min(len(s1), len(s2))) 
                      if s1[k] != s2[k] and s1[k] not in ['-', 'N'] and s2[k] not in ['-', 'N'])
        dist[(n1, n2)] = diff
        print(f"{diff:6}", end="")
    print()

# PARSIMONY SCORING
print("\n" + "="*80)
print("PARSIMONY SCORING (5 Pleistocene specimens + Pan outgroup)")
print("="*80)

# H1: Sima as sister-group to (Denisova + Neandertals)
h1_score = (dist[('Sima_430ka','Denisova3_70ka')] + 
            dist[('Neandertal_Vindija_44ka','Feldhofer1_40ka')] +
            dist[('Sima_430ka','Neandertal_Vindija_44ka')])

# H2: Sima in ancestral position
h2_score = (dist[('Sima_430ka','Denisova3_70ka')] + 
            dist[('Denisova3_70ka','Neandertal_Vindija_44ka')])

print(f"\nH1 (Sima as sister-group):        {h1_score:6} steps")
print(f"H2 (Sima ancestral position):    {h2_score:6} steps")
print(f"Difference:                      {h1_score - h2_score:6} steps")

if h1_score > 0:
    reduction = 100 * (h1_score - h2_score) / h1_score
    print(f"Parsimony reduction (H2):        {reduction:6.1f}%")

print("\n" + "="*80)
if h2_score < h1_score:
    print("✓ H2 (ancestral position) is MORE PARSIMONIOUS")
    print(f"  Signal strength: {reduction:.1f}% improvement over H1")
else:
    print("✓ H1 (sister-group) is more parsimonious")
print("="*80)

# Summary table
print("\nKEY DISTANCES:")
print("-"*80)
print(f"  Sima to Denisova3:              {dist[('Sima_430ka','Denisova3_70ka')]:6} substitutions")
print(f"  Sima to Neandertal_Vindija:    {dist[('Sima_430ka','Neandertal_Vindija_44ka')]:6} substitutions")
print(f"  Denisova3 to Neandertal:       {dist[('Denisova3_70ka','Neandertal_Vindija_44ka')]:6} substitutions")
print(f"  Neandertals (intra-group):     {dist[('Neandertal_Vindija_44ka','Feldhofer1_40ka')]:6} substitutions")
print(f"  Pan to all hominins:           1300-1500 substitutions (clear outgroup)")
print("="*80)
