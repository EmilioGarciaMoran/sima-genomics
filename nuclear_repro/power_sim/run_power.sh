#!/usr/bin/env bash
# Simulación de potencia completa: 3 escenarios de T_ND, 100 réplicas por condición.
#   micromamba activate powersim && bash run_power.sh [hilos]
# Tiempo aproximado con 2 núcleos: ~25 min. Las simulaciones (sims/, ~100 MB) no se versionan.
set -euo pipefail
cd "$(dirname "$0")"
T=${1:-2}
command -v python3 >/dev/null && python3 -c "import msprime, numpy" 2>/dev/null || {
  echo "ERROR: activa el entorno: micromamba activate powersim" >&2; exit 1; }

declare -A SEED=( [470]=470001 [440]=440001 [550]=550001 )
for TND in 470 440 550; do echo "$TND ${SEED[$TND]}"; done | \
  xargs -P "$T" -n 2 sh -c 'python3 scripts/01_simulate.py "$0" "$1" "sims/TND$0" > "sims_TND$0.log" 2>&1 && echo "simulado TND$0"'

for TND in 470 440 550; do
  python3 scripts/02_power.py "sims/TND$TND" "TND$TND" "results/TND$TND" 100
done
python3 scripts/03_summary.py results/TND*/D_REPLICATES.tsv | tee results/SUMMARY.tsv
for TND in 470 440 550; do cp "sims/TND$TND/meta.json" "results/TND$TND/sim_meta.json"; done
rm -f sims_TND*.log
