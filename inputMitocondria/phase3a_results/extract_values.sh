#!/bin/bash

echo "================================================================================"
echo "PHASE 3a: FINAL RESULTS EXTRACTION"
echo "================================================================================"
echo ""

echo "IQ-TREE (HKY+G4) - Maximum Likelihood"
echo "======================================"
cat iqtree_hky.iqtree | grep "Log-likelihood" | head -5
echo ""
cat iqtree_hky.treefile
echo ""

echo "RAxML-NG (GTR+G) - Maximum Likelihood"
echo "====================================="
cat raxml_gtr.raxml.log | grep "Final GAMMA"
echo ""
cat raxml_gtr.raxml.support
echo ""

echo "================================================================================"

