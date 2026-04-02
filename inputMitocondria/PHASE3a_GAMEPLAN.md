# PHASE 3a GAME PLAN: ML/BAYES OPTIMIZATION OF mtDNA PHYLOGENY
## Local execution on egarmo-Dell-3431 with IQ-TREE, RAxML-NG, MrBayes

**Date:** March 21, 2026  
**Machine:** egarmo-Dell-3431 (i7-9700, 8 cores, 15GB RAM)  
**Software versions:** IQ-TREE 3.0.1, RAxML-NG 1.2.2, MrBayes 3.2.7

---

## 📋 WORKFLOW OVERVIEW

```
INPUT: alignment_6magnificos.fasta (6 mtDNA genomes, 16,568 bp)
  ↓
[STEP 1] DATA PREP: FASTA → PHYLIP (15 min)
  ↓
[STEP 2] IQ-TREE ML + Bootstrap (5 min)
  ↓
[STEP 3] RAxML-NG ML + Bootstrap (10 min)
  ↓
[STEP 4] MrBayes MCMC (30-60 min)
  ↓
[STEP 5] KH Topology Tests (5-10 min)
  ↓
OUTPUT: Support values + Tree topologies + Statistical tests
  ↓
[STEP 6] Parse results → Update AJBA manuscript
```

---

## 🚀 EXECUTION INSTRUCTIONS

### Setup Phase

```bash
# 1. Download the three scripts to your local machine:
#    - phase3a_master.sh
#    - parse_results.py
#    - compare_topologies.py

# 2. Copy scripts to working directory
cd ~/sima/inputMitocondria
cp /path/to/phase3a_master.sh .
cp /path/to/parse_results.py .
cp /path/to/compare_topologies.py .

# 3. Make scripts executable
chmod +x phase3a_master.sh
chmod +x parse_results.py
chmod +x compare_topologies.py

# 4. Activate phylo environment
micromamba activate phylo

# 5. Verify all software is installed
which iqtree2 raxml-ng mb
```

### Analysis Phase

```bash
# 6. Run the master analysis workflow
# This will take 45-90 minutes depending on your CPU
bash phase3a_master.sh

# Output will be in: ~/sima/inputMitocondria/phase3a_results/
```

### Results Extraction Phase

```bash
# 7. Parse and summarize results
cd phase3a_results
python3 ../parse_results.py

# 8. Compare topologies across methods
python3 ../compare_topologies.py

# This generates summary tables showing:
# - Bootstrap support values (IQ-TREE, RAxML)
# - Posterior probabilities (MrBayes)
# - KH/AU topology test p-values
# - Which topology is better supported
```

---

## 📊 EXPECTED OUTPUTS

### IQ-TREE (HKY+G4)
- `iqtree_hky.treefile` - Best ML tree with bootstrap values
- `iqtree_hky.iqtree` - Full analysis log with statistics
- Expected: Bootstrap ≥95% for well-supported clades

### RAxML-NG (GTR+G)
- `raxml_gtr.raxml.bestTree` - Best ML tree
- `raxml_gtr.raxml.support` - Tree with bootstrap support values
- Expected: Bootstrap ≥95% for well-supported clades

### MrBayes (HKY+G)
- `mrbayes.con.tre` - 50% majority rule consensus tree
- `mrbayes.p` - Posterior samples (chain file)
- `mrbayes.sum` - Summary statistics
- Expected: Posterior probability ≥0.95 for well-supported clades

### KH Topology Tests
- `khtest_h0_mtdna.iqtree` - Test results for mtDNA topology (H0)
- `khtest_h1_ndna.iqtree` - Test results for nDNA topology (H1)
- AU test p-value interpretation:
  - p > 0.05: Topology NOT significantly rejected
  - p ≤ 0.05: Topology significantly rejected

---

## 🔍 KEY RESULTS TO EXTRACT

Once analysis completes, extract these values for the AJBA manuscript:

```
1. BOOTSTRAP SUPPORT (IQ-TREE & RAxML)
   - Value for Sima basal position
   - Value for (Denisova + Neandertal) clade
   - Value for ((Feldhofer + Mezmaiskaya) clade

2. POSTERIOR PROBABILITY (MrBayes)
   - Same clades as above

3. TOPOLOGY TESTS (KH)
   - AU p-value for H0 (mtDNA topology)
   - AU p-value for H1 (nDNA topology)
   - Which hypothesis is better supported?

4. TREE LIKELIHOOD
   - Log-likelihood from IQ-TREE
   - Log-likelihood from RAxML
   - Log-likelihood from KH tests
```

---

## 📝 MANUSCRIPT INTEGRATION

### New Methods section entries:

```
Maximum likelihood phylogenetic reconstruction was performed using 
IQ-TREE v3.0.1 (model HKY+G4) and RAxML-NG v1.2.2 (model GTR+G), 
each with 1,000 bootstrap replicates. Bayesian MCMC inference was 
conducted in MrBayes v3.2.7 with 1,000,000 generations (HKY+G model). 
Topology testing was performed using the Kishino-Hasegawa test 
(10,000 replicates) as implemented in IQ-TREE, evaluating two competing 
hypotheses: (H0) Sima basal to (Denisova + Neandertal), and (H1) 
(Sima + Neandertal) sister to Denisova.
```

### New Results section entries:

```
Both IQ-TREE (bootstrap [VALUE]%) and RAxML-NG (bootstrap [VALUE]%) 
strongly support the basal position of Sima mtDNA relative to the 
(Denisova + Neandertal) clade. MrBayes posterior probability for this 
topology is [VALUE]. The Kishino-Hasegawa test significantly rejects 
the alternative (nDNA) topology (H1; AU p = [VALUE]), providing 
independent statistical support for the mtDNA-derived topology.
```

---

## ⏱️ TIMING ESTIMATES

| Step | Time | Notes |
|------|------|-------|
| Data prep | 2-5 min | FASTA→PHYLIP conversion |
| IQ-TREE | 5 min | 1000 bootstrap, HKY+G4 |
| RAxML-NG | 10 min | 1000 bootstrap, GTR+G |
| MrBayes | 45-60 min | 1M generations, 8 cores |
| KH tests | 5-10 min | 10,000 replicates each test |
| **TOTAL** | **65-95 min** | Can run overnight |

---

## 🛠️ TROUBLESHOOTING

### IQ-TREE not found
```bash
micromamba install -c bioconda iqtree -y
```

### RAxML-NG not found
```bash
micromamba install -c bioconda raxml-ng -y
```

### MrBayes not found
```bash
micromamba install -c bioconda mrbayes -y
```

### PHYLIP format error
- Check that sequence names are ≤10 characters
- All sequences must be same length
- Run: `python3 parse_results.py` to validate

### MrBayes hangs
- Ctrl+C to stop
- Reduce `ngen` in the mbcommands.txt to 100,000 for testing
- Or run overnight with default settings

---

## 📋 CHECKLIST

- [ ] Scripts downloaded and placed in ~/sima/inputMitocondria
- [ ] Scripts made executable (chmod +x)
- [ ] phylo environment activated
- [ ] IQ-TREE, RAxML-NG, MrBayes verified (which command)
- [ ] alignment_6magnificos.fasta exists in working directory
- [ ] phase3a_master.sh execution started
- [ ] Results appear in phase3a_results/ subdirectory
- [ ] parse_results.py executed successfully
- [ ] compare_topologies.py executed successfully
- [ ] Key values extracted and noted for manuscript
- [ ] Methods section updated
- [ ] Results section updated

---

## 🔗 NEXT STEPS IN CLAUDE

Once you have run all scripts and extracted results:

1. **Post in Claude chat:**
   - Directory listing of phase3a_results/
   - Output from parse_results.py
   - Output from compare_topologies.py

2. **I will:**
   - Verify all analyses completed correctly
   - Extract exact numbers for the manuscript
   - Calculate Cohen's d effect sizes if needed
   - Rewrite Methods section with exact parameters
   - Rewrite Results section with exact support values
   - Update Discussion with ML/Bayesian findings

3. **Final output:**
   - Updated AJBA manuscript with mtDNA ML/Bayes section complete
   - Figures showing tree topologies with support values
   - Tables comparing methods

---

## 📚 REFERENCES FOR OUTPUT INTERPRETATION

- **IQ-TREE:** https://doi.org/10.1093/molbev/msy148
- **RAxML-NG:** https://doi.org/10.1093/bioinformatics/btz305
- **MrBayes:** https://doi.org/10.1093/bioinformatics/btg180
- **KH test:** https://doi.org/10.1093/oxfordjournals.molbev.a025679

---

**Questions?** Ping me in Claude after you run the analysis!

🚀 Ready to execute Phase 3a!
