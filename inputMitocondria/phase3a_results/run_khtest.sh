#!/bin/bash
# KH Test Script - Ejecutar en: ~/sima/inputMitocondria/phase3a_results

cd ~/sima/inputMitocondria/phase3a_results

# 1. Crear constraint trees con nombres EXACTOS
cat > constraint_h0_mtdna.tre << 'EOF'
(Sima_430ka,(Denisova3_70ka,(Neandertal_Vindija_44ka,(Feldhofer1_40ka,Mezmaiskaya1_60ka))));
EOF

cat > constraint_h1_ndna.tre << 'EOF'
((Sima_430ka,Denisova3_70ka),(Neandertal_Vindija_44ka,(Feldhofer1_40ka,Mezmaiskaya1_60ka)));
EOF

echo "✅ Constraint trees created"

# 2. Activar conda
micromamba activate phylo

# 3. Ejecutar KH tests
echo "🔄 Running KH test H0 (basal Sima)..."
iqtree -s alignment_6magnificos.phy -m HKY+G4 -z constraint_h0_mtdna.tre -zb 10000 -au -nt AUTO -pre khtest_h0_mtdna 2>&1 | tee khtest_h0_mtdna_run.log

echo "🔄 Running KH test H1 (sister Sima-Denisova)..."
iqtree -s alignment_6magnificos.phy -m HKY+G4 -z constraint_h1_ndna.tre -zb 10000 -au -nt AUTO -pre khtest_h1_ndna 2>&1 | tee khtest_h1_ndna_run.log

# 4. Extraer p-values
echo ""
echo "========================================"
echo "RESULTADOS KH TEST"
echo "========================================"
echo ""
echo "H0 (basal Sima):"
grep -A 5 "USER TREES\|p-value\|AU test" khtest_h0_mtdna.iqtree | head -20

echo ""
echo "H1 (sister Sima-Denisova):"
grep -A 5 "USER TREES\|p-value\|AU test" khtest_h1_ndna.iqtree | head -20
