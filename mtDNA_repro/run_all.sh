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

# Comprobar el entorno antes de empezar (evita, p. ej., que 'python' sea Python 2 del sistema)
missing=""
for c in mafft iqtree3 raxml-ng python3; do command -v "$c" >/dev/null || missing="$missing $c"; done
python3 -c "import sys, Bio, matplotlib; assert sys.version_info >= (3, 8)" 2>/dev/null || missing="$missing biopython/matplotlib(python3)"
if [ -n "$missing" ]; then
  echo "ERROR: faltan herramientas:$missing" >&2
  echo "Activa el entorno:  micromamba activate mtphylo   (o créalo: micromamba create -f environment.yml)" >&2
  exit 1
fi
bash scripts/01_get_sequences.sh "$MODE"
bash scripts/02_align.sh
bash scripts/03_trees.sh "$T"
bash scripts/04_topology_tests.sh "$T"
bash scripts/05_legacy_check.sh
python3 scripts/06_figure.py es
python3 scripts/06_figure.py en
echo
echo "Hecho. Compara sequences/SEQUENCES.tsv y los results/*.tsv con la versión del repositorio (git diff)."
