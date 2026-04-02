#!/usr/bin/env python3
"""
Phase 2: Compile Master FASTA (6 Magnificent mtDNA)
Input: Individual FASTA files from phase2_download_CORRECTED.py
Output: pleistocene_6magnificos_MASTER_CORRECTED.fasta
"""
import os
from Bio import SeqIO

os.chdir('inputMitocondria')

print("="*80)
print("PHASE 2: COMPILE MASTER FASTA (CORRECTED)")
print("="*80)

# Expected files (from corrected download)
specimens = {
    'Sima_430ka': 'Sima_femurXIII_430ka_KF683087.fasta',
    'Denisova3_70ka': 'Denisova3_70ka_NC_013993.fasta',
    'Neandertal_Vindija_44ka': 'Neandertal_Vindija_44ka_NC_020627.fasta',
    'Feldhofer1_40ka': 'Feldhofer1_40ka_FM865407.fasta',
    'Mezmaiskaya1_60ka': 'Mezmaiskaya1_60ka_FM865411.fasta',
    'Pan_outgroup': 'Pan_troglodytes_NC_001643.fasta',
}

sequences = {}
found = 0
missing = 0

print("\nCompiling sequences...\n")

for label, filename in specimens.items():
    if os.path.exists(filename):
        try:
            for record in SeqIO.parse(filename, 'fasta'):
                sequences[label] = str(record.seq).upper()
                print(f"  ✓ {label:30} {len(record.seq):6} bp")
                found += 1
        except Exception as e:
            print(f"  ✗ {label:30} ERROR: {str(e)[:40]}")
            missing += 1
    else:
        print(f"  ✗ {label:30} FILE NOT FOUND: {filename}")
        missing += 1

# Write master FASTA
output_file = 'pleistocene_6magnificos_MASTER_CORRECTED.fasta'

if found == len(specimens):
    with open(output_file, 'w') as f:
        for label, seq in sequences.items():
            f.write(f">{label}\n{seq}\n")
    
    total_bp = sum(len(s) for s in sequences.values())
    print("\n" + "="*80)
    print(f"✓ MASTER FASTA CREATED: {output_file}")
    print(f"  Sequences: {len(sequences)}")
    print(f"  Total bp: {total_bp}")
    print("="*80)
    print("\nNext: run MAFFT alignment")
    print(f"  mafft --thread 4 {output_file} > alignment_6magnificos_CORRECTED.fasta")
    
else:
    print("\n" + "="*80)
    print(f"⚠️  INCOMPLETE: {found}/{len(specimens)} sequences found")
    print("="*80)
