#!/usr/bin/env bash
# Paso 4: test de topologías. Las tres hipótesis se evalúan en UNA misma ejecución,
# así IQ-TREE puede calcular p-valores comparativos (KH, SH, AU). En la fase 3a cada
# hipótesis se evaluó en ejecuciones separadas, con un solo árbol cada una, y por eso
# nunca hubo p-valores.
#
#   H0: Sima basal a (Denisova + Neandertales)   <- hipótesis del preprint
#   H1: Sima hermano de Denisova                  <- Meyer et al. 2014
#   H2: Sima hermano de los Neandertales
#
#   bash scripts/04_topology_tests.sh [hilos]
#
# Modelos: HKY+G4 (el del preprint) y el elegido por BIC en el paso 3, sobre full y coding.
# Salidas: results/topology/* y results/TOPOLOGY_TESTS.tsv
# Ejecutar desde mtDNA_repro/ después del paso 3.
set -euo pipefail
T=${1:-2}; [ "$T" -gt "$(nproc)" ] && T=$(nproc)   # IQ-TREE aborta si T > núcleos
SEED=12345
OUT=results/topology
mkdir -p "$OUT"

cat > "$OUT/hypotheses.trees" <<'EOF'
(Pan_outgro,(Sima_430ka,(Denisova3_,(Mezmaiskay,(Neandertal,Feldhofer1)))));
(Pan_outgro,((Sima_430ka,Denisova3_),(Mezmaiskay,(Neandertal,Feldhofer1))));
(Pan_outgro,(Denisova3_,(Sima_430ka,(Mezmaiskay,(Neandertal,Feldhofer1)))));
EOF

S=results/TOPOLOGY_TESTS.tsv
printf "alineamiento\tmodelo\thipotesis\tlogL\tdeltaL\tp_KH\tp_SH\tp_AU\n" > "$S"
for A in full coding; do
  BIC=$(grep -m1 "Best-fit model according to BIC" "results/trees/iq_${A}_mfp_ufboot.iqtree" | awk '{print $NF}')
  for M in HKY+G4 "$BIC"; do
    P="$OUT/${A}_${M//+/_}"
    iqtree3 -s "alignments/aln_${A}.fasta" -m "$M" -o Pan_outgro -z "$OUT/hypotheses.trees" \
            -zb 10000 -au -seed $SEED -T 1 -pre "$P" -redo -quiet   # 1 hilo: p-valores idénticos en cualquier máquina
    sed -n '/^Tree      logL/,/^$/p' "$P.iqtree" | awk -v a="$A" -v m="$M" '
      $1 ~ /^[123]$/ {
        h = ($1==1) ? "H0_Sima_basal" : ($1==2) ? "H1_Sima+Denisova" : "H2_Sima+Neandertales"
        # columnas: n logL deltaL bp-RELL [+-] p-KH [+-] p-SH [+-] c-ELW [+-] p-AU [+-]
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", a, m, h, $2, $3, $6, $8, $12
      }' >> "$S"
  done
done
cat "$S"
