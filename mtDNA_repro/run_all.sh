#!/usr/bin/env bash
# Ejecuta el pipeline mtDNA completo, de las secuencias a los tests de topología.
#
#   micromamba create -f environment.yml && micromamba activate mtphylo
#   bash run_all.sh [genbank|recover] [hilos]
#
# genbank (por defecto): descarga las secuencias de NCBI. recover: las extrae del alineamiento legacy.
# Duración aproximada con 2 núcleos: ~16 min (el bootstrap estándar de IQ-TREE es lo más lento).
set -euo pipefail
cd "$(dirname "$0")"
MODE=${1:-genbank}
T=${2:-2}
bash scripts/01_get_sequences.sh "$MODE"
bash scripts/02_align.sh
bash scripts/03_trees.sh "$T"
bash scripts/04_topology_tests.sh "$T"
bash scripts/05_legacy_check.sh
echo
echo "Hecho. Compara sequences/SEQUENCES.tsv y los results/*.tsv con la versión del repositorio (git diff)."
