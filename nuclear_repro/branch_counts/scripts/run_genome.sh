#!/usr/bin/env bash
# Recorrido genómico: para cada cromosoma, descarga los cinco VCF, extrae los sitios
# informativos y borra los VCF antes de pasar al siguiente. Pico de disco ~6 GB.
#
#   bash scripts/run_genome.sh [dir_trabajo] [primer_cr] [ultimo_cr]
#   nohup bash scripts/run_genome.sh > results/genoma.log 2>&1 &
#
# REANUDABLE: si un cromosoma ya tiene su .npz en results/genome/, se salta. Si se corta
# a mitad, vuelve a lanzarlo y sigue donde iba.
#
# Salida: results/genome/chrN_inf.npz, uno por cromosoma. Son pequeños (decenas de MB):
# solo contienen posiciones con Mbuti ancestral y al menos un arcaico derivado, que es el
# superconjunto de todos los conjuntos de referencia.
set -uo pipefail
cd "$(dirname "$0")/.."

WORK=${1:-/tmp/sima_vcf}
FIRST=${2:-1}
LAST=${3:-22}
OUT=results/genome
mkdir -p "$OUT" "$WORK"

echo "=================================================================="
echo "Recorrido genómico  |  inicio: $(date '+%F %T')"
echo "  cromosomas $FIRST a $LAST"
echo "  VCF temporales en: $WORK (se borran tras cada cromosoma)"
echo "  salidas en: $OUT"
echo "=================================================================="

libre() { df -BG --output=avail "$WORK" 2>/dev/null | tail -1 | tr -dc '0-9'; }

ok=0; fallos=""
for C in $(seq "$FIRST" "$LAST"); do
  NPZ="$OUT/chr${C}_inf.npz"
  if [ -s "$NPZ" ]; then
    echo "[cr$C] ya hecho, se salta"
    ok=$((ok+1)); continue
  fi

  L=$(libre)
  if [ -n "$L" ] && [ "$L" -lt 12 ]; then
    echo "[cr$C] ABORTA: solo ${L} GB libres en $WORK" >&2
    fallos="$fallos $C"; break
  fi

  echo
  echo "------------------------------------------------------------------"
  echo "[cr$C] $(date '+%T')  descargando..."
  rm -f "$WORK"/*_chr*.vcf.gz "$WORK"/*_chr*.vcf.gz.tbi
  if ! bash scripts/fetch_chrom.sh "$C" "$WORK"; then
    echo "[cr$C] FALLO en la descarga" >&2
    fallos="$fallos $C"; continue
  fi

  echo "[cr$C] $(date '+%T')  procesando..."
  if python3 scripts/read_vcfs.py --dir "$WORK" --chrom "$C" \
       --out "$NPZ" --mode homozygous --only-informative --progress 0; then
    ls -lh "$NPZ" | awk -v c="$C" '{print "  cr"c": "$5" en "$9}'
    ok=$((ok+1))
  else
    echo "[cr$C] FALLO al procesar" >&2
    rm -f "$NPZ"
    fallos="$fallos $C"
  fi

  rm -f "$WORK"/*_chr*.vcf.gz "$WORK"/*_chr*.vcf.gz.tbi
  echo "[cr$C] $(date '+%T')  hecho. VCF borrados."
done

echo
echo "=================================================================="
echo "Fin: $(date '+%F %T')   cromosomas completados: $ok"
[ -n "$fallos" ] && echo "FALLARON O QUEDAN PENDIENTES:$fallos  (relanza el script: es reanudable)"
du -sh "$OUT" 2>/dev/null
echo "=================================================================="
