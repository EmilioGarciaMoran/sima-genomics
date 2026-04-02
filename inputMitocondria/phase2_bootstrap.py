#!/usr/bin/env python3
"""
PHASE 2: BOOTSTRAP VALIDATION
Test robustness of H1 topology (Sima sister-group)

Run from: ~/sima/inputMitocondria/
Input: alignment_6magnificos.fasta

Usage: python3 phase2_bootstrap.py
"""
import random
from Bio import SeqIO
import os

os.chdir('.')

print("="*80)
print("PHASE 2: BOOTSTRAP VALIDATION (1000 replicates)")
print("="*80)

alignment_file = 'alignment_6magnificos.fasta'

if not os.path.exists(alignment_file):
    print(f"ERROR: {alignment_file} not found")
    exit(1)

# Read alignment
print(f"\nReading alignment: {alignment_file}")
seqs = {}
for record in SeqIO.parse(alignment_file, 'fasta'):
    seqs[record.id] = str(record.seq).upper()

align_length = len(list(seqs.values())[0])
print(f"  Sequences: {len(seqs)}")
print(f"  Alignment length: {align_length} bp")

# Define taxa
taxa = [
    ('Sima_430ka', seqs.get('Sima_430ka')),
    ('Denisova3_70ka', seqs.get('Denisova3_70ka')),
    ('Neandertal_Vindija_44ka', seqs.get('Neandertal_Vindija_44ka')),
    ('Feldhofer1_40ka', seqs.get('Feldhofer1_40ka')),
    ('Mezmaiskaya1_60ka', seqs.get('Mezmaiskaya1_60ka')),
    ('Pan_outgroup', seqs.get('Pan_outgroup')),
]

# Verify all present
missing = [name for name, seq in taxa if seq is None]
if missing:
    print(f"ERROR: Missing sequences: {missing}")
    exit(1)

# Bootstrap function
def compute_parsimony(taxa_list, alignment_positions):
    """Compute parsimony scores for H1 and H2 using given positions"""
    
    # Build subsequences from selected positions
    sub_taxa = []
    for name, seq in taxa_list:
        subseq = ''.join(seq[i] for i in alignment_positions if i < len(seq))
        sub_taxa.append((name, subseq))
    
    # Pairwise distances
    dist = {}
    for i, (n1, s1) in enumerate(sub_taxa):
        for j, (n2, s2) in enumerate(sub_taxa):
            if i == j:
                diff = 0
            else:
                diff = sum(1 for k in range(min(len(s1), len(s2))) 
                          if s1[k] != s2[k] and s1[k] not in ['-', 'N'] and s2[k] not in ['-', 'N'])
            dist[(n1, n2)] = diff
    
    # Parsimony H1: Sima sister-group
    h1_score = (dist[('Sima_430ka','Denisova3_70ka')] + 
                dist[('Neandertal_Vindija_44ka','Feldhofer1_40ka')] +
                dist[('Sima_430ka','Neandertal_Vindija_44ka')])
    
    # Parsimony H2: Sima ancestral
    h2_score = (dist[('Sima_430ka','Denisova3_70ka')] + 
                dist[('Denisova3_70ka','Neandertal_Vindija_44ka')])
    
    return h1_score, h2_score

# Run bootstrap
n_bootstrap = 1000
h1_wins = 0
h2_wins = 0
h1_scores = []
h2_scores = []

print(f"\nRunning {n_bootstrap} bootstrap replicates...")
print("-"*80)

for rep in range(n_bootstrap):
    # Resample positions with replacement
    positions = [random.randint(0, align_length - 1) for _ in range(align_length)]
    
    # Compute parsimony on resampled alignment
    h1, h2 = compute_parsimony(taxa, positions)
    
    h1_scores.append(h1)
    h2_scores.append(h2)
    
    if h1 < h2:
        h1_wins += 1
    else:
        h2_wins += 1
    
    # Progress
    if (rep + 1) % 100 == 0:
        print(f"  {rep + 1:4d}/{n_bootstrap} replicates...")

print("\n" + "="*80)
print("BOOTSTRAP RESULTS")
print("="*80)

h1_pct = 100 * h1_wins / n_bootstrap
h2_pct = 100 * h2_wins / n_bootstrap

print(f"\nH1 (Sima sister-group):  {h1_wins:4d}/{n_bootstrap}  ({h1_pct:5.1f}%)")
print(f"H2 (Sima ancestral):     {h2_wins:4d}/{n_bootstrap}  ({h2_pct:5.1f}%)")

# Statistics
import statistics

h1_mean = statistics.mean(h1_scores)
h1_stdev = statistics.stdev(h1_scores)
h2_mean = statistics.mean(h2_scores)
h2_stdev = statistics.stdev(h2_scores)

print(f"\nH1 scores: mean={h1_mean:.1f} ± {h1_stdev:.1f}")
print(f"H2 scores: mean={h2_mean:.1f} ± {h2_stdev:.1f}")

print("\n" + "="*80)
print("INTERPRETATION")
print("="*80)

if h1_pct >= 95:
    print(f"\n✓ H1 is STRONGLY SUPPORTED ({h1_pct:.1f}% bootstrap support)")
    print("  High confidence in sister-group topology")
elif h1_pct >= 75:
    print(f"\n✓ H1 is MODERATELY SUPPORTED ({h1_pct:.1f}% bootstrap support)")
    print("  Reasonable confidence but some uncertainty")
elif h1_pct >= 50:
    print(f"\n⚠ H1 is WEAKLY SUPPORTED ({h1_pct:.1f}% bootstrap support)")
    print("  Weak signal; consider alternative explanations")
else:
    print(f"\n✗ H2 is actually PREFERRED ({h2_pct:.1f}% bootstrap support)")
    print("  H1 may be artifact of sampling")

print("\n" + "="*80)
