#!/usr/bin/env python3
"""
PHASE 2: CLEAN MTDNA ANALYSIS
6 Magnificent - Verified accessions only
Run from: ~/sima/inputMitocondria/

Usage: python3 phase2_clean.py
"""
import subprocess
import os
from Bio import SeqIO

print("="*80)
print("PHASE 2: CLEAN MTDNA ANALYSIS (6 MAGNIFICENT)")
print("="*80)
print()

# ACCESSIONS (verified by DeepSeek + ChatGPT)
accessions = {
    'Sima_430ka': 'KF683087',
    'Denisova3_70ka': 'NC_013993',         # ✓ VERIFIED mtDNA
    'Neandertal_Vindija_44ka': 'AM948965', # ✓ VERIFIED mtDNA
    'Feldhofer1_40ka': 'FM865407',
    'Mezmaiskaya1_60ka': 'FM865411',
    'Pan_outgroup': 'NC_001643',
}

# ============================================================================
# STEP 1: DOWNLOAD
# ============================================================================
print("STEP 1: DOWNLOADING SEQUENCES")
print("-"*80)

downloaded = {}

for name, acc in accessions.items():
    print(f"  {name:30} ({acc:12})...", end=' ', flush=True)
    try:
        result = subprocess.run(
            ['/usr/bin/efetch', '-db', 'nucleotide', '-id', acc, '-format', 'fasta'],
            capture_output=True, text=True, timeout=30
        )
        
        if result.returncode == 0 and len(result.stdout) > 100:
            filename = f"{name}_{acc}.fasta"
            with open(filename, 'w') as f:
                f.write(result.stdout)
            
            lines = result.stdout.strip().split('\n')
            seq_lines = [l for l in lines[1:] if not l.startswith('>')]
            seq_len = sum(len(l) for l in seq_lines)
            
            print(f"✓ {seq_len:6} bp")
            downloaded[name] = (acc, seq_len, filename)
        else:
            print(f"✗ FAILED")
    except Exception as e:
        print(f"✗ ERROR")

if len(downloaded) < len(accessions):
    print(f"\n✗ ERROR: Only {len(downloaded)}/{len(accessions)} downloaded")
    exit(1)

print(f"\n✓ Downloaded: {len(downloaded)}/{len(accessions)}")

# ============================================================================
# STEP 2: COMPILE MASTER FASTA
# ============================================================================
print("\nSTEP 2: COMPILE MASTER FASTA")
print("-"*80)

sequences = {}
for name, (acc, bp, filename) in downloaded.items():
    for record in SeqIO.parse(filename, 'fasta'):
        sequences[name] = str(record.seq).upper()
        print(f"  ✓ {name:30} {bp:6} bp")

master_file = 'pleistocene_6magnificos_master.fasta'
with open(master_file, 'w') as f:
    for label, seq in sequences.items():
        f.write(f">{label}\n{seq}\n")

total_bp = sum(len(s) for s in sequences.values())
print(f"\n✓ Master FASTA: {master_file} ({total_bp} bp total)")

# ============================================================================
# STEP 3: MAFFT ALIGNMENT
# ============================================================================
print("\nSTEP 3: MAFFT ALIGNMENT")
print("-"*80)

align_file = 'alignment_6magnificos.fasta'
print(f"  Input: {master_file}")
print(f"  Output: {align_file}")
print(f"  Running: mafft --thread 4...\n")

try:
    with open(align_file, 'w') as out:
        result = subprocess.run(
            ['mafft', '--thread', '4', master_file],
            stdout=out, stderr=subprocess.PIPE, text=True, timeout=300
        )
    
    if result.returncode == 0:
        lines = open(align_file).readlines()
        n_seqs = sum(1 for line in lines if line.startswith('>'))
        print(f"✓ Alignment complete: {len(lines)} lines, {n_seqs} sequences")
    else:
        print(f"✗ MAFFT failed")
        exit(1)
except Exception as e:
    print(f"✗ ERROR: {str(e)}")
    exit(1)

# ============================================================================
# STEP 4: PARSIMONY ANALYSIS
# ============================================================================
print("\nSTEP 4: PARSIMONY ANALYSIS")
print("-"*80)

def read_alignment(filename):
    seqs = {}
    for record in SeqIO.parse(filename, 'fasta'):
        seqs[record.id] = str(record.seq).upper()
    return seqs

seqs = read_alignment(align_file)

taxa = [
    ('Sima_430ka', seqs.get('Sima_430ka')),
    ('Denisova3_70ka', seqs.get('Denisova3_70ka')),
    ('Neandertal_Vindija_44ka', seqs.get('Neandertal_Vindija_44ka')),
    ('Feldhofer1_40ka', seqs.get('Feldhofer1_40ka')),
    ('Mezmaiskaya1_60ka', seqs.get('Mezmaiskaya1_60ka')),
    ('Pan_outgroup', seqs.get('Pan_outgroup')),
]

# Verify all sequences
missing = [name for name, seq in taxa if seq is None]
if missing:
    print(f"✗ ERROR: Missing sequences: {missing}")
    exit(1)

print("\nPAIRWISE DISTANCE MATRIX:")
print("-"*80)

dist = {}
for i, (n1, s1) in enumerate(taxa):
    print(f"{n1:25}", end="")
    for j, (n2, s2) in enumerate(taxa):
        if i == j:
            diff = 0
        else:
            diff = sum(1 for k in range(min(len(s1), len(s2))) 
                      if s1[k] != s2[k] and s1[k] not in ['-', 'N'] and s2[k] not in ['-', 'N'])
        dist[(n1, n2)] = diff
        print(f"{diff:6}", end="")
    print()

# Parsimony
h1_score = (dist[('Sima_430ka','Denisova3_70ka')] + 
            dist[('Neandertal_Vindija_44ka','Feldhofer1_40ka')] +
            dist[('Sima_430ka','Neandertal_Vindija_44ka')])

h2_score = (dist[('Sima_430ka','Denisova3_70ka')] + 
            dist[('Denisova3_70ka','Neandertal_Vindija_44ka')])

print("\n" + "="*80)
print("PARSIMONY RESULTS")
print("="*80)

print(f"\nH1 (Sima sister-group):          {h1_score:6} steps")
print(f"H2 (Sima ancestral):             {h2_score:6} steps")
print(f"Difference:                      {h1_score - h2_score:6} steps")

if h1_score > 0:
    reduction = 100 * (h1_score - h2_score) / h1_score
    print(f"Reduction (H2):                  {reduction:6.1f}%")

print("\n" + "="*80)
if h2_score < h1_score:
    print("✓ RESULT: H2 (ancestral) is MORE PARSIMONIOUS")
else:
    print("✓ RESULT: H1 (sister-group) is more parsimonious")
print("="*80)

print("\nOUTPUT FILES:")
print(f"  {master_file}")
print(f"  {align_file}")
print("\n✓ PHASE 2 COMPLETE")
