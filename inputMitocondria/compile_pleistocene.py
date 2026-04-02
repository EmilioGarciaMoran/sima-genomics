#!/usr/bin/env python3
from Bio import SeqIO
import os

os.chdir('inputMitocondria')

sequences = {}

files_to_add = [
    'Sima_KF683087.fasta',
    'Denisova2_mtdna.fasta',
    'Denisova3_FN673705.fasta',
    'Denisova8_mtdna.fasta',
    'Denisova15_mtdna.fasta',
    'Denisova19_mtdna.fasta',
    'Denisova20_mtdna.fasta',
    'Denisova21_mtdna.fasta',
    'Neandertal_mtdna.fasta',
    'HomoSapiens_mtdna.fasta'
]

for fasta_file in files_to_add:
    if os.path.exists(fasta_file):
        for record in SeqIO.parse(fasta_file, 'fasta'):
            name = record.id.split('|')[0]
            sequences[name] = str(record.seq).upper()
            print(f"✓ {name}: {len(record.seq)} bp")
    else:
        print(f"⚠ {fasta_file} not found")

with open('pleistocene_master.fasta', 'w') as f:
    for name, seq in sequences.items():
        f.write(f">{name}\n{seq}\n")

print(f"\n✓ Master FASTA compiled: {len(sequences)} sequences")
print(f"  Total bp: {sum(len(s) for s in sequences.values())}")

