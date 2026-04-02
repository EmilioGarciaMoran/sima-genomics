#!/bin/bash
# Phase 2: Download ALL Pleistocene mtDNA from GenBank
# Core MSA: Sima + Denisova (2,3,8,15,19,20,21) + Neandertal + Homo sapiens
# External: Ust-Ishim, Pan (pairwise post-hoc)

cd inputMitocondria

echo "======================================"
echo "DOWNLOADING COMPLETE PLEISTOCENE SET"
echo "======================================"

# DENISOVA
echo "Denisova 2 (~195 ka)..."
efetch -db nucleotide -id "NC_024829" -format fasta > Denisova2_mtdna.fasta

echo "Denisova 3 (FN673705, already have)..."
# Ya tienes

echo "Denisova 8 (KT780370)..."
efetch -db nucleotide -id "KT780370" -format fasta > Denisova8_mtdna.fasta

echo "Denisova 15 (MK033602)..."
efetch -db nucleotide -id "MK033602" -format fasta > Denisova15_mtdna.fasta

echo "Denisova 19 (MT576650)..."
efetch -db nucleotide -id "MT576650" -format fasta > Denisova19_mtdna.fasta

echo "Denisova 20 (MT576651)..."
efetch -db nucleotide -id "MT576651" -format fasta > Denisova20_mtdna.fasta

echo "Denisova 21 (MT576652)..."
efetch -db nucleotide -id "MT576652" -format fasta > Denisova21_mtdna.fasta

# NEANDERTAL
echo "Neandertal (EU260658)..."
efetch -db nucleotide -id "EU260658" -format fasta > Neandertal_mtdna.fasta

# HOMO SAPIENS
echo "Homo sapiens (NC_012920)..."
efetch -db nucleotide -id "NC_012920" -format fasta > HomoSapiens_mtdna.fasta

# OUTGROUP/EXTERNAL
echo "Pan troglodytes (NC_001643, outgroup)..."
efetch -db nucleotide -id "NC_001643" -format fasta > Pan_mtdna.fasta

echo "======================================"
echo "✓ All downloads complete!"
ls -lh *_mtdna.fasta
echo "======================================"
