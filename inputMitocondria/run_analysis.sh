#!/bin/bash
set -e

cd phase3a_results

echo "====== IQ-TREE ======"
iqtree -s alignment_6magnificos.phy -m HKY+G4 -bb 1000 -nt AUTO -seed 12345 -pre iqtree_hky -redo

echo ""
echo "====== RAxML-NG ======"
raxml-ng --msa alignment_6magnificos.phy --model GTR+G --threads AUTO --seed 12345 --prefix raxml_gtr --all --bs-trees 1000

echo ""
echo "====== KH Test H0 ======"
iqtree -s alignment_6magnificos.phy -m HKY+G4 -z constraint_h0_mtdna.tre -zb 10000 -au -nt AUTO -pre khtest_h0_mtdna

echo ""
echo "====== KH Test H1 ======"
iqtree -s alignment_6magnificos.phy -m HKY+G4 -z constraint_h1_ndna.tre -zb 10000 -au -nt AUTO -pre khtest_h1_ndna

echo ""
echo "✅ DONE"

