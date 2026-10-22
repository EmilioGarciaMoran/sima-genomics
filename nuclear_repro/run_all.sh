#!/usr/bin/env bash
# Ejecuta el pipeline nuclear sobre los BAM recortados (slices/), sin AWS ni genoma completo.
#
#   micromamba create -f environment.yml && micromamba activate nucrepro
#   bash run_all.sh                        # datos reales versionados en slices/
#   bash run_all.sh SLICES CATALOG RESULTS # otros directorios (los tests usan datos sintéticos)
#
# Los BAM recortados se generan una sola vez, donde están los BAM completos:
#   python3 scripts/00_make_slices.py --ref hg19.fa --bam /ruta/*.bam
set -euo pipefail
cd "$(dirname "$0")"
SLICES=${1:-slices}
CATALOG=${2:-catalog}
RESULTS=${3:-results}

python3 -c "import sys, pysam; assert sys.version_info >= (3, 9)" 2>/dev/null || {
  echo "ERROR: hace falta python3 >= 3.9 con pysam. Activa el entorno: micromamba activate nucrepro" >&2
  exit 1
}

echo "== 1. Comprobación de entradas";            python3 scripts/01_check_inputs.py "$SLICES" "$CATALOG"
echo "== 2. Genotipado directo";                  python3 scripts/02_direct_genotype.py "$SLICES" "$CATALOG" "$RESULTS"
echo "== 3. Verosimilitudes genotípicas (CDS)";   python3 scripts/03_genotype_likelihoods.py "$SLICES" "$CATALOG" "$RESULTS"
echo "== 3b. ANGSD (opcional)";                   bash scripts/03b_angsd.sh "$SLICES" "$CATALOG" "$RESULTS"
echo "== 4. Informe";                             python3 scripts/04_report.py "$SLICES" "$RESULTS"
echo
echo "Hecho. Compara $RESULTS/*.tsv con la versión del repositorio (git diff)."
