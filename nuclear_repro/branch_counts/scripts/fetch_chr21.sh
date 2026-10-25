#!/usr/bin/env bash
# Descarga los VCF del cromosoma 21 necesarios para la validación (ver DATOS_NECESARIOS.md).
# Unos 2,5 GB en total. Reanudable: wget -c salta lo ya descargado.
#
#   bash scripts/fetch_chr21.sh [directorio_destino]
#
# No descarga el estado ancestral de Ensembl: comprobar antes el nombre exacto del fichero.
set -euo pipefail
DEST=${1:-data/chr21}
BASE=https://cdna.eva.mpg.de/neandertal
mkdir -p "$DEST"
cd "$DEST"

command -v wget >/dev/null || { echo "ERROR: hace falta wget" >&2; exit 1; }

get() {  # get <subruta> <nombre_local>
  echo "--> $2"
  wget -c -q --show-progress -O "$2" "$BASE/$1" || { echo "FALLO en $1" >&2; return 1; }
}

get "Vindija/VCF/Altai/chr21_mq25_mapab100.vcf.gz"          "D5_chr21.vcf.gz"
get "Vindija/VCF/Altai/chr21_mq25_mapab100.vcf.gz.tbi"      "D5_chr21.vcf.gz.tbi"
get "Vindija/VCF/Vindija33.19/chr21_mq25_mapab100.vcf.gz"      "Vi3319_chr21.vcf.gz"
get "Vindija/VCF/Vindija33.19/chr21_mq25_mapab100.vcf.gz.tbi"  "Vi3319_chr21.vcf.gz.tbi"
get "Vindija/VCF/Denisova/chr21_mq25_mapab100.vcf.gz"       "D3_chr21.vcf.gz"
get "Vindija/VCF/Denisova/chr21_mq25_mapab100.vcf.gz.tbi"   "D3_chr21.vcf.gz.tbi"
get "Chagyrskaya/VCF/chr21.noRB.vcf.gz"                     "Chag8_chr21.vcf.gz"
get "Chagyrskaya/VCF/chr21.noRB.vcf.gz.tbi"                 "Chag8_chr21.vcf.gz.tbi"
get "altai/ModernHumans/vcf/SS6004471.hg19_1000g.21.mod.vcf.gz"     "Mbuti_chr21.vcf.gz"
get "altai/ModernHumans/vcf/SS6004471.hg19_1000g.21.mod.vcf.gz.tbi" "Mbuti_chr21.vcf.gz.tbi"

# md5 de origen, para comprobar a mano (los nombres locales difieren de los remotos)
for d in Vindija/VCF/Altai Vindija/VCF/Vindija33.19 Vindija/VCF/Denisova Chagyrskaya/VCF; do
  wget -q -O "md5_$(echo "$d" | tr / _).txt" "$BASE/$d/md5sums.txt" 2>/dev/null || true
done

echo
echo "Descargado en $(pwd):"
ls -lh *.vcf.gz
echo
echo "Comprueba los md5 contra los ficheros md5_*.txt (el nombre remoto es chr21_*.vcf.gz"
echo "o chr21.noRB.vcf.gz, el local lleva el nombre del espécimen)."
