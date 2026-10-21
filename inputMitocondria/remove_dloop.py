#!/usr/bin/env python3
# DEPRECATED (2026-10-20): NO USAR.
# Este script elimina las primeras 1122 columnas del alineamiento, que no son el D-loop.
# En coordenadas rCRS el D-loop ocupa 16024-16569 y 1-576 (final y principio de la molécula).
# Sustituido por mtDNA_repro/scripts/02_align.sh (aln_coding.fasta).

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

def write_fasta(filename, seqs):
    with open(filename, 'w') as f:
        for name, seq in seqs.items():
            f.write('>' + name + '\n' + seq + '\n')

print("mtDNA ANALYSIS - CODING REGIONS ONLY (NO D-LOOP)")
print("="*60)

seqs = read_fasta('alignment_6_mtdna.fasta')

dloop_size = 1122

print("\nOriginal alignment stats:")
for name, seq in seqs.items():
    print("  " + name[:20] + " " + str(len(seq)) + " bp")
    break

print("\nRemoving D-loop region (first " + str(dloop_size) + " bp)")

seqs_coding = {}
for name, seq in seqs.items():
    coding_seq = seq[dloop_size:]
    seqs_coding[name] = coding_seq

print("\nCoding regions stats:")
for name, seq in seqs_coding.items():
    print("  " + name[:20] + " " + str(len(seq)) + " bp")
    break

write_fasta('alignment_6_mtdna_CODING.fasta', seqs_coding)
print("\nSaved: alignment_6_mtdna_CODING.fasta")

sima = None
d3 = None
nean = None
hsap = None
pan = None

for k, v in seqs_coding.items():
    if 'KF683087' in k:
        sima = v
