#!/usr/bin/env bash
# Paso 3: árboles de máxima verosimilitud y bootstrap, enraizados EXPLÍCITAMENTE en Pan.
#
#   bash scripts/03_trees.sh [hilos]
#
# Para cada alineamiento (full y coding):
#   a) IQ-TREE: selección de modelo (BIC) + árbol ML + UFBoot 1000 (réplicas guardadas)
#   b) IQ-TREE: HKY+G4 (el modelo del preprint) + bootstrap estándar no paramétrico 1000
#   c) RAxML-NG: GTR+G + bootstrap estándar 1000
# Después cuenta en las réplicas el soporte de cada agrupación (clade_support.py).
# Salidas: results/trees/* y results/CLADE_SUPPORT.tsv
# Ejecutar desde mtDNA_repro/ con el entorno mtphylo activo.
set -euo pipefail
T=${1:-2}; [ "$T" -gt "$(nproc)" ] && T=$(nproc)   # IQ-TREE aborta si T > núcleos
SEED=12345
OUT=results/trees
mkdir -p "$OUT"
iqtree3 --version | awk 'NR==1' > "$OUT/VERSIONS.txt"
raxml-ng --version | awk '/RAxML/ && !f {print; f=1}' >> "$OUT/VERSIONS.txt"

for A in full coding; do
  ALN=alignments/aln_${A}.fasta
  # UFBoot con un solo hilo: con varios hilos el resultado cambia ligeramente según su número
  # (comprobado: 4 frente a 2 hilos dieron 97,8 frente a 98,0). Con -T 1 es igual en cualquier máquina.
  iqtree3 -s "$ALN" -m MFP -o Pan_outgro -B 1000 -wbt -seed $SEED -T 1 -pre "$OUT/iq_${A}_mfp_ufboot" -redo -quiet
  iqtree3 -s "$ALN" -m HKY+G4 -o Pan_outgro -b 1000 -seed $SEED -T "$T" -pre "$OUT/iq_${A}_hky_stdboot" -redo -quiet
  raxml-ng --all --msa "$ALN" --model GTR+G --outgroup Pan_outgro --bs-trees 1000 \
           --seed $SEED --threads "$T" --prefix "$OUT/rx_${A}_gtr" --redo > /dev/null
done

S=results/CLADE_SUPPORT.tsv
python3 -I scripts/clade_support.py --header > "$S"
for A in full coding; do
  M=$(grep -m1 "Best-fit model according to BIC" "$OUT/iq_${A}_mfp_ufboot.iqtree" | awk '{print $NF}')
  python3 -I scripts/clade_support.py "$OUT/iq_${A}_mfp_ufboot.ufboot"   "${A}_IQTREE_${M}_UFBoot" >> "$S"
  python3 -I scripts/clade_support.py "$OUT/iq_${A}_hky_stdboot.boottrees" "${A}_IQTREE_HKY+G4_stdboot" >> "$S"
  python3 -I scripts/clade_support.py "$OUT/rx_${A}_gtr.raxml.bootstraps" "${A}_RAxML_GTR+G_stdboot" >> "$S"
done

# Árbol ML enraizado de cada análisis, para inspección
for f in "$OUT"/iq_*_mfp_ufboot.contree "$OUT"/iq_*_hky_stdboot.contree "$OUT"/rx_*_gtr.raxml.support; do
  printf "%s\t%s\n" "$(basename "$f")" "$(cat "$f")"
done > results/ML_TREES.tsv

cat "$S"
