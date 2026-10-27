#!/usr/bin/env bash
# Pieza 3: extrae los alelos de Sima en los sitios de referencia del genoma completo.
#
#   bash scripts/extract_sima.sh [directorio_datos]
#
# Reanudable: salta el espécimen cuyo npz ya exista.
set -uo pipefail
cd "$(dirname "$0")/.."

DATA=${1:-$HOME/sima/sima_data}
BAM="$DATA/bam_dmg"
SITIOS=../branch_counts/results/genoma_inf.npz
OUT=results/alelos
mkdir -p "$OUT"

[ -s "$SITIOS" ] || { echo "ERROR: falta $SITIOS (ejecuta merge_genome.py)" >&2; exit 1; }
command -v samtools >/dev/null || { echo "ERROR: falta samtools" >&2; exit 1; }

echo "=================================================================="
echo "Extracción de alelos  |  inicio: $(date '+%F %T')"
echo "  sitios: $SITIOS"
echo "=================================================================="

saca() {  # saca <etiqueta> <bam>
  local etq=$1 bam=$2
  local sal="$OUT/${etq}.npz"
  [ -s "$sal" ] && { echo "  [$etq] ya hecho"; return 0; }
  [ -s "$bam" ] || { echo "  [$etq] falta $bam" >&2; return 1; }
  python3 scripts/extract_alleles.py --sitios "$SITIOS" --bam "$bam" \
          --etiqueta "$etq" --out "$sal" || { rm -f "$sal"; return 1; }
}

while IFS=$'\t' read -r esp muestra corridas lecturas bases papel; do
  case "$esp" in \#*|especimen|"") continue ;; esac
  echo
  echo "------------------------------------------------------------------"
  echo "[$esp] $(date '+%T')  ($papel)"
  saca "$esp" "$BAM/${esp}.bam"
  # variantes de sensibilidad, solo para los dos principales y el de confirmación
  case "$papel" in
    principal|confirmacion)
      for v in nofilt ss ds; do
        [ -s "$BAM/${esp}.${v}.bam" ] && saca "${esp}.${v}" "$BAM/${esp}.${v}.bam"
      done ;;
  esac
done < ESPECIMENES.tsv

echo
echo "=================================================================="
echo "Fin: $(date '+%F %T')"
printf "  %-26s %12s %12s\n" fichero "sitios" derivados
for f in "$OUT"/*.npz; do
  [ -s "$f" ] || continue
  python3 - "$f" <<'PY'
import sys, numpy as np, os
z = np.load(sys.argv[1]); c = z["codigo"]
print(f"  {os.path.basename(sys.argv[1])[:-4]:<26s} "
      f"{int((c>=0).sum()):12,} {int((c==1).sum()):12,}")
PY
done
echo
echo "Siguiente y última pieza: recuento por rama y D-statistics."
echo "=================================================================="
