#!/bin/bash
################################################################################
# PHASE 3a: ML/BAYES OPTIMIZATION OF mtDNA PHYLOGENY
# Local execution on egarmo-Dell-3431 with IQ-TREE, RAxML-NG, MrBayes
# 
# Usage: bash phase3a_master.sh
# 
# Workflow:
# 1. DATA PREP: FASTA → PHYLIP
# 2. ML ANALYSIS: IQ-TREE + Bootstrap
# 3. ML ANALYSIS: RAxML-NG + Bootstrap
# 4. BAYESIAN: MrBayes MCMC
# 5. TOPOLOGY TEST: KH test (IQ-TREE)
################################################################################

set -e

# CONFIGURATION
WORKDIR="$(pwd)"
FASTA_INPUT="alignment_6magnificos.fasta"
PHYLIP_INPUT="alignment_6magnificos.phy"
OUTPUT_DIR="phase3a_results"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")

echo "================================================================================"
echo "PHASE 3a: ML/BAYES mtDNA PHYLOGENY"
echo "Start time: $(date)"
echo "Working directory: $WORKDIR"
echo "================================================================================"
echo ""

# CREATE OUTPUT DIRECTORY
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

################################################################################
# STEP 1: DATA PREPARATION (FASTA → PHYLIP)
################################################################################
echo "[STEP 1] DATA PREPARATION: FASTA → PHYLIP"
echo "=========================================="
echo ""

# Check input
if [ ! -f "../$FASTA_INPUT" ]; then
    echo "❌ ERROR: $FASTA_INPUT not found in parent directory"
    exit 1
fi

echo "Converting $FASTA_INPUT to PHYLIP format..."

# Python script inline for conversion
python3 << 'PYTHON_EOF'
# Read FASTA
fasta_file = "../alignment_6magnificos.fasta"
records = []
current_id = None
current_seq = []

with open(fasta_file) as f:
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            if current_id:
                records.append((current_id, ''.join(current_seq)))
            current_id = line[1:]
            current_seq = []
        else:
            current_seq.append(line)
    if current_id:
        records.append((current_id, ''.join(current_seq)))

n_seqs = len(records)
seq_len = len(records[0][1])

print(f"Sequences: {n_seqs}")
print(f"Length: {seq_len} bp")
print()

# Write PHYLIP (sequential format)
with open("alignment_6magnificos.phy", "w") as f:
    f.write(f"{n_seqs} {seq_len}\n")
    for seq_id, seq in records:
        # Truncate/pad name to 10 chars (PHYLIP standard)
        name = seq_id[:10].ljust(10)
        f.write(f"{name}{seq}\n")

print(f"✅ PHYLIP file created: alignment_6magnificos.phy")
PYTHON_EOF

echo ""

################################################################################
# STEP 2: IQ-TREE ML + BOOTSTRAP
################################################################################
echo "[STEP 2] IQ-TREE MAXIMUM LIKELIHOOD + 1000 BOOTSTRAP"
echo "====================================================="
echo ""
echo "Model: HKY+G4"
echo "Bootstrap replicates: 1000"
echo "Threads: AUTO"
echo ""

iqtree2 \
    -s "$PHYLIP_INPUT" \
    -m HKY+G4 \
    -bb 1000 \
    -nt AUTO \
    -seed 12345 \
    -pre iqtree_hky \
    -redo \
    2>&1 | tee iqtree_hky.log

echo ""
echo "✅ IQ-TREE completed"
echo "   Best tree: iqtree_hky.treefile"
echo "   Bootstrap: 1000 replicates with support values"
echo ""

################################################################################
# STEP 3: RAxML-NG ML + BOOTSTRAP
################################################################################
echo "[STEP 3] RAxML-NG MAXIMUM LIKELIHOOD + 1000 BOOTSTRAP"
echo "======================================================"
echo ""
echo "Model: GTR+G"
echo "Bootstrap replicates: 1000"
echo ""

raxml-ng \
    --msa "$PHYLIP_INPUT" \
    --model GTR+G \
    --threads AUTO \
    --seed 12345 \
    --prefix raxml_gtr \
    --all \
    --bs-trees 1000 \
    2>&1 | tee raxml_gtr.log

echo ""
echo "✅ RAxML-NG completed"
echo "   Best tree: raxml_gtr.raxml.bestTree"
echo "   Bootstrap: 1000 replicates"
echo ""

################################################################################
# STEP 4: MrBayes BAYESIAN INFERENCE (MCMC)
################################################################################
echo "[STEP 4] MrBayes BAYESIAN INFERENCE (MCMC)"
echo "=========================================="
echo ""
echo "Model: HKY+G"
echo "MCMC chains: 1,000,000 generations"
echo "Burnin: 250,000 generations (25%)"
echo "Sampling frequency: 100"
echo ""

# Create MrBayes command file
cat > mbcommands.txt << 'MRBAYES_EOF'
set autoclose=yes nowarn=yes;
execute alignment_6magnificos.phy;
lset nst=2 rates=gamma;
mcmc ngen=1000000 samplefreq=100 printfreq=1000 diagnfreq=5000 filename=mrbayes;
sump burnin=2500;
sumt burnin=2500 contype=halfcompat;
quit;
MRBAYES_EOF

echo "Running MrBayes..."
mb < mbcommands.txt 2>&1 | tee mrbayes.log

echo ""
echo "✅ MrBayes MCMC completed"
echo "   Posterior samples: mrbayes.p"
echo "   Consensus tree (50%): mrbayes.con.tre"
echo ""

################################################################################
# STEP 5: KISHINO-HASEGAWA (KH) TOPOLOGY TEST
################################################################################
echo "[STEP 5] KISHINO-HASEGAWA (KH) TOPOLOGY TEST"
echo "============================================="
echo ""
echo "Comparing two competing topologies:"
echo "  H0 (mtDNA): Sima basal sister to (Denisova + Neandertal)"
echo "  H1 (nDNA):  (Sima + Neandertal) sister to Denisova"
echo ""

# Create constraint tree for H0 (mtDNA topology)
cat > constraint_h0_mtdna.tre << 'TREE_EOF'
(Sima_430ka,(Denisova3_70ka,(Neandertal_Vindija_44ka,(Feldhofer1_40ka,Mezmaiskaya1_60ka))));
TREE_EOF

# Create constraint tree for H1 (nDNA topology)
cat > constraint_h1_ndna.tre << 'TREE_EOF'
((Sima_430ka,Neandertal_Vindija_44ka),(Denisova3_70ka,(Feldhofer1_40ka,Mezmaiskaya1_60ka)));
TREE_EOF

echo "Testing H0 (mtDNA topology)..."
iqtree2 \
    -s "$PHYLIP_INPUT" \
    -m HKY+G4 \
    -z constraint_h0_mtdna.tre \
    -zb 10000 \
    -au \
    -nt AUTO \
    -pre khtest_h0_mtdna \
    2>&1 | tee khtest_h0_mtdna.log

echo ""
echo "Testing H1 (nDNA topology)..."
iqtree2 \
    -s "$PHYLIP_INPUT" \
    -m HKY+G4 \
    -z constraint_h1_ndna.tre \
    -zb 10000 \
    -au \
    -nt AUTO \
    -pre khtest_h1_ndna \
    2>&1 | tee khtest_h1_ndna.log

echo ""
echo "✅ KH Tests completed"
echo "   H0 results: khtest_h0_mtdna.iqtree"
echo "   H1 results: khtest_h1_ndna.iqtree"
echo ""

################################################################################
# STEP 6: SUMMARY & RESULTS EXTRACTION
################################################################################
echo "[STEP 6] EXTRACTING KEY RESULTS"
echo "==============================="
echo ""

python3 << 'RESULTS_EOF'
import re

print("=" * 80)
print("PHASE 3a RESULTS SUMMARY")
print("=" * 80)
print()

# Extract IQ-TREE bootstrap values
print("IQ-TREE (HKY+G4) Results:")
print("-" * 40)
try:
    with open("iqtree_hky.iqtree") as f:
        content = f.read()
        # Look for log likelihood
        if "Log-likelihood of the tree:" in content:
            for line in content.split('\n'):
                if "Log-likelihood" in line:
                    print(f"  {line.strip()}")
        # Count bootstrap replicates
        if "1000 bootstrap trees" in content:
            print("  ✅ Bootstrap: 1000 replicates completed")
except:
    print("  (waiting for results)")

print()

# Extract RAxML results
print("RAxML-NG (GTR+G) Results:")
print("-" * 40)
try:
    with open("raxml_gtr.raxml.info") as f:
        for line in f:
            if "best-fit model" in line.lower() or "log likelihood" in line.lower():
                print(f"  {line.strip()}")
except:
    print("  (waiting for results)")

print()

# Extract MrBayes summary
print("MrBayes (HKY+G) Results:")
print("-" * 40)
try:
    with open("mrbayes.sum") as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if "Clade" in line or "Probability" in line:
                print(f"  {line.strip()}")
except:
    print("  (waiting for results)")

print()
print("=" * 80)
print("NEXT STEP: Run parse_results.py to extract bootstrap/posterior values")
print("=" * 80)

RESULTS_EOF

echo ""

################################################################################
# FINAL SUMMARY
################################################################################
echo ""
echo "================================================================================"
echo "✅ PHASE 3a WORKFLOW COMPLETED"
echo "================================================================================"
echo ""
echo "Output directory: $(pwd)"
echo ""
echo "Files generated:"
echo "  IQ-TREE:     iqtree_hky.treefile, iqtree_hky.iqtree"
echo "  RAxML-NG:    raxml_gtr.raxml.bestTree, raxml_gtr.raxml.support"
echo "  MrBayes:     mrbayes.con.tre, mrbayes.p"
echo "  KH tests:    khtest_h0_mtdna.iqtree, khtest_h1_ndna.iqtree"
echo ""
echo "NEXT: Parse results with parse_results.py and update AJBA manuscript"
echo ""
