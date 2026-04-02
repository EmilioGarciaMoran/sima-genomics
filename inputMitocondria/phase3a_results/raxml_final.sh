#!/bin/bash
################################################################################
# RAXML-NG FINAL EXECUTION
# Convierte PHYLIP a FASTA y ejecuta RAxML-NG con bootstrap
#
# Usage: bash raxml_final.sh
# Run from: ~/sima/inputMitocondria/phase3a_results/
################################################################################

set -e

echo "================================================================================"
echo "RAXML-NG: MAXIMUM LIKELIHOOD + 1000 BOOTSTRAP"
echo "Start time: $(date)"
echo "Working directory: $(pwd)"
echo "================================================================================"
echo ""

# Step 1: Convert PHYLIP to FASTA
echo "[STEP 1] Converting PHYLIP to FASTA format..."
echo "=========================================="
echo ""

python3 << 'PYEOF'
print("Reading PHYLIP file...")
with open("alignment_6magnificos.phy") as f:
    lines = f.readlines()

header = lines[0].split()
n_seqs = int(header[0])
seq_len = int(header[1])

print(f"Sequences: {n_seqs}, Length: {seq_len} bp")

sequences = {}
for line in lines[1:]:
    if line.strip():
        parts = line.split()
        name = parts[0]
        seq = parts[1]
        sequences[name] = seq

print(f"Parsed {len(sequences)} sequences")
print()

print("Writing FASTA file...")
with open("alignment_6magnificos.fasta", "w") as f:
    for name, seq in sequences.items():
        f.write(f">{name}\n{seq}\n")

print("✅ FASTA file created: alignment_6magnificos.fasta")
print(f"   Sequences: {len(sequences)}")
print(f"   Length: {len(list(sequences.values())[0])} bp")

PYEOF

echo ""

# Step 2: Run RAxML-NG
echo "[STEP 2] Running RAxML-NG (GTR+G model, 1000 bootstrap replicates)..."
echo "=========================================="
echo ""

raxml-ng \
    --msa alignment_6magnificos.fasta \
    --msa-format FASTA \
    --model GTR+G \
    --threads AUTO \
    --seed 12345 \
    --prefix raxml_gtr \
    --all \
    --bs-trees 1000 \
    2>&1 | tee raxml_gtr_final.log

echo ""
echo "================================================================================"
echo "✅ RAxML-NG ANALYSIS COMPLETED"
echo "================================================================================"
echo ""
echo "Output files:"
echo "  Best ML tree:              raxml_gtr.raxml.bestTree"
echo "  Tree with bootstrap:       raxml_gtr.raxml.support"
echo "  Full analysis log:         raxml_gtr.raxml.log"
echo "  Session log:               raxml_gtr_final.log"
echo ""
echo "Next: Review results and compare topologies"
echo ""
