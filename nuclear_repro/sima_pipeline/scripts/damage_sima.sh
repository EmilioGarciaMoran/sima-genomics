#!/usr/bin/env bash
# Pieza 2: filtro de daño y enmascarado, aplicado a los BAM de la pieza 1.
#
#   bash scripts/damage_sima.sh [directorio_datos] [hilos]
#
# Produce <especimen>.dmg.bam, y además las variantes de los análisis de sensibilidad
# preregistrados (apartado 10 ter): sin filtro de daño, y convención de doble cadena.
# Reanudable: salta el que ya exista.
set -uo pipefail
cd "$(dirname "$0")/.."

DATA=${1:-$HOME/sima/sima_data}
T=${2:-4}
REF="$DATA/ref/human_g1k_v37.fasta"
BAM="$DATA/bam"
OUT="$DATA/bam_dmg"
mkdir -p "$OUT" results

echo "=================================================================="
echo "Filtro de daño  |  inicio: $(date '+%F %T')"
echo "  criterio: C->T o G->A en la primera o ultima posicion alineada (Meyer 2016)"
echo "  enmascarado: todas las timinas en las 3 primeras y 3 ultimas posiciones"
echo "  salida: $OUT"
echo "=================================================================="

for t in samtools; do
  command -v "$t" >/dev/null || { echo "ERROR: falta $t" >&2; exit 1; }
done
[ -s "$REF.fai" ] || { echo "ERROR: falta $REF.fai (samtools faidx)" >&2; exit 1; }

# variante <sufijo> <argumentos extra>
corre() {
  local esp=$1 suf=$2; shift 2
  local sal="$OUT/${esp}${suf}.bam"
  [ -s "$sal" ] && { echo "  [$esp$suf] ya hecho"; return 0; }
  samtools view -h "$BAM/${esp}.bam" \
    | python3 scripts/damage_filter.py --ref "$REF" --etiqueta "$esp$suf" "$@" \
    | samtools sort -@ "$T" -o "$sal" - 2>/dev/null \
    && samtools index "$sal" \
    || { echo "  [$esp$suf] FALLO" >&2; rm -f "$sal"; return 1; }
}

while IFS=$'\t' read -r esp muestra corridas lecturas bases papel; do
  case "$esp" in \#*|especimen|"") continue ;; esac
  [ -s "$BAM/${esp}.bam" ] || { echo "[$esp] falta el BAM de la pieza 1" >&2; continue; }
  echo
  echo "------------------------------------------------------------------"
  echo "[$esp] $(date '+%T')  ($papel)  Meyer declara $bases bases útiles"
  corre "$esp" ""        --library any --mask 3 --mask-mode T
  corre "$esp" ".nofilt" --library any --mask 3 --mask-mode T --no-filter
  corre "$esp" ".ss"     --library ss  --mask 3 --mask-mode hebra
  corre "$esp" ".ds"     --library ds  --mask 3 --mask-mode hebra
done < ESPECIMENES.tsv

echo
echo "=================================================================="
echo "Fin: $(date '+%F %T')"
echo
printf "  %-20s %12s %12s %10s %10s\n" especimen conservadas "tasa daño" "Meyer" "cociente"
while IFS=$'\t' read -r esp muestra corridas lecturas bases papel; do
  case "$esp" in \#*|especimen|"") continue ;; esac
  b="$OUT/${esp}.bam"; [ -s "$b" ] || continue
  ent=$(samtools view -c "$BAM/${esp}.bam")
  sal=$(samtools view -c "$b")
  tasa=$(awk -v a="$sal" -v b="$ent" 'BEGIN{if(b>0)printf "%.1f%%",100*a/b; else print "-"}')
  printf "  %-20s %12s %12s %10s %10s\n" "$esp" "$sal" "$tasa" "$bases" "-"
done < ESPECIMENES.tsv
echo
echo "La tasa de daño debe ser ALTA en los especímenes endógenos y BAJA en la escápula,"
echo "contaminada al 98%: el ADN moderno no lleva desaminación terminal."
echo "=================================================================="
