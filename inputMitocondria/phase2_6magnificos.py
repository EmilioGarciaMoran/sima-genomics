#!/usr/bin/env python3
"""
Phase 2: The 6 Magnificent - Complete Pleistocene mtDNA MSA
Verified, complete genomes only (~16.5 kb each)
"""
import subprocess
import os
from Bio import SeqIO

os.chdir('inputMitocondria')

accessions = {
    'Sima_femurXIII': 'KF683087',
    'Denisova3': 'FN673705',
    'Neandertal_Altai': 'KC879692',
    'Feldhofer1': 'FM865407',
    'Mezmaiskaya1': 'FM865411',
    'Pan_troglodytes': 'NC_001643',
}

print("="*70)
print("PHASE 2: THE 6 MAGNIFICENT - Complete mtDNA")
print("="*70)

downloaded = {}

for name, accession in accessions.items():
    print(f"  {name:25} ({accession})...", end=' ', flush=True)
    result = subprocess.run(
        ['/usr/bin/efetch', '-db', 'nucleotide', '-id', accession, '-format', 'fasta'],
        capture_output=True, text=True, timeout=30
    )
    
    if result.returncode == 0 and len(result.stdout) > 100:
        filename = f"{name}_{accession}.fasta"
        with open(filename, 'w') as f:
            f.write(result.stdout)
        
        lines = result.stdout.strip().split('\n')
        seq_lines = [l for l in lines[1:] if not l.startswith('>')]
        seq_len = sum(len(l) for l in seq_lines)
        print(f"✓ {seq_len} bp")
        downloaded[name] = seq_len

print("\n" + "="*70)
print(f"Downloaded: {len(downloaded)}/6")
print("="*70)

# Compile master FASTA
print("\nCompiling master FASTA...")
sequences = {}
for name, accession in accessions.items():
    for record in SeqIO.parse(f"{name}_{accession}.fasta", 'fasta'):
        sequences[name] = str(record.seq).upper()

with open('pleistocene_6magnificos_master.fasta', 'w') as f:
    for name, seq in sequences.items():
        f.write(f">{name}\n{seq}\n")

print(f"✓ Master FASTA: {len(sequences)} sequences")
print(f"  Total bp: {sum(len(s) for s in sequences.values())}")
