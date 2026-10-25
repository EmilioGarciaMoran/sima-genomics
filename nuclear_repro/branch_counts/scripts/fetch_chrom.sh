#!/usr/bin/env bash
# Descarga los cinco VCF de UN cromosoma. Generaliza fetch_chr21.sh.
#
#   bash scripts/fetch_chrom.sh <cromosoma> <directorio_destino>
#
# Comprueba las 10 URL antes de descargar: si falta alguna, no baja nada y sale con error.
set -euo pipefail

C=${1:?uso: fetch_chrom.sh <cromosoma> <destino>}
DEST=${2:?uso: fetch_chrom.sh <cromosoma> <destino>}
HOSTS=(https://ftp.eva.mpg.de/neandertal https://cdna.eva.mpg.de/neandertal)

FILES=$(cat <<EOF
D5_chr$C.vcf.gz	Vindija/VCF/Altai/chr${C}_mq25_mapab100.vcf.gz
D5_chr$C.vcf.gz.tbi	Vindija/VCF/Altai/chr${C}_mq25_mapab100.vcf.gz.tbi
Vi3319_chr$C.vcf.gz	Vindija/VCF/Vindija33.19/chr${C}_mq25_mapab100.vcf.gz
Vi3319_chr$C.vcf.gz.tbi	Vindija/VCF/Vindija33.19/chr${C}_mq25_mapab100.vcf.gz.tbi
D3_chr$C.vcf.gz	Vindija/VCF/Denisova/chr${C}_mq25_mapab100.vcf.gz
D3_chr$C.vcf.gz.tbi	Vindija/VCF/Denisova/chr${C}_mq25_mapab100.vcf.gz.tbi
Chag8_chr$C.vcf.gz	Chagyrskaya/VCF/chr${C}.noRB.vcf.gz
Chag8_chr$C.vcf.gz.tbi	Chagyrskaya/VCF/chr${C}.noRB.vcf.gz.tbi
Mbuti_chr$C.vcf.gz	altai/ModernHumans/vcf/SS6004471.hg19_1000g.${C}.mod.vcf.gz
Mbuti_chr$C.vcf.gz.tbi	altai/ModernHumans/vcf/SS6004471.hg19_1000g.${C}.mod.vcf.gz.tbi
EOF
)

BASE=""
for h in "${HOSTS[@]}"; do
  curl -sIf --max-time 25 "$h/Vindija/VCF/Altai/chr${C}_mq25_mapab100.vcf.gz" >/dev/null 2>&1 \
    && { BASE=$h; break; }
done
[ -n "$BASE" ] || { echo "ERROR cr$C: ningún servidor responde" >&2; exit 1; }

missing=0
while IFS=$'\t' read -r local remote; do
  [ -n "$local" ] || continue
  curl -sIf --max-time 25 "$BASE/$remote" >/dev/null 2>&1 || {
    echo "  FALTA: $remote" >&2; missing=$((missing+1)); }
done <<< "$FILES"
[ "$missing" -eq 0 ] || { echo "ERROR cr$C: $missing fichero(s) no encontrados" >&2; exit 1; }

mkdir -p "$DEST"; cd "$DEST"
while IFS=$'\t' read -r local remote; do
  [ -n "$local" ] || continue
  want=$(curl -sIf --max-time 25 "$BASE/$remote" \
         | awk 'BEGIN{IGNORECASE=1}/^content-length:/{gsub(/\r/,"");print $2}' | tail -1)
  wget -c -q -O "$local" "$BASE/$remote" || { echo "ERROR cr$C: descarga de $remote" >&2; exit 1; }
  got=$(stat -c%s "$local" 2>/dev/null || echo 0)
  # COMPROBACIÓN IMPRESCINDIBLE: un .gz truncado no da error al descomprimir con zcat,
  # el lector lo procesa como si estuviera completo y escribe resultados parciales
  # que parecen válidos. Esto ya ocurrió una vez por falta de espacio en disco.
  if [ -n "$want" ] && [ "$got" != "$want" ]; then
    echo "ERROR cr$C: $local incompleto ($got de $want bytes)" >&2
    rm -f "$local"; exit 1
  fi
  case "$local" in
    *.vcf.gz) gzip -t "$local" 2>/dev/null || {
        echo "ERROR cr$C: $local no supera la prueba de integridad gzip" >&2
        rm -f "$local"; exit 1; } ;;
  esac
done <<< "$FILES"
du -sh . | awk -v c="$C" '{print "  cr"c": descargado "$1" (tamaños e integridad verificados)"}'
