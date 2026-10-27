#!/usr/bin/env bash
# Descarga los VCF del cromosoma 21 necesarios para la validación (ver DATOS_NECESARIOS.md).
# Unos 2,5 GB. Reanudable: vuelve a lanzarlo y continúa donde se quedó.
#
#   bash scripts/fetch_chr21.sh [directorio_destino]
#
# Primero comprueba que las 10 URL existen (HEAD). Si alguna falla, no descarga nada:
# así un cambio en el servidor se detecta antes de bajar gigabytes.
#
# No descarga el estado ancestral de Ensembl: ver DATOS_NECESARIOS.md.
set -euo pipefail

DEST=${1:-data/chr21}
HOSTS=(https://ftp.eva.mpg.de/neandertal https://cdna.eva.mpg.de/neandertal)

# nombre_local <TAB> subruta remota
FILES=$(cat <<'EOF'
D5_chr21.vcf.gz	Vindija/VCF/Altai/chr21_mq25_mapab100.vcf.gz
D5_chr21.vcf.gz.tbi	Vindija/VCF/Altai/chr21_mq25_mapab100.vcf.gz.tbi
Vi3319_chr21.vcf.gz	Vindija/VCF/Vindija33.19/chr21_mq25_mapab100.vcf.gz
Vi3319_chr21.vcf.gz.tbi	Vindija/VCF/Vindija33.19/chr21_mq25_mapab100.vcf.gz.tbi
D3_chr21.vcf.gz	Vindija/VCF/Denisova/chr21_mq25_mapab100.vcf.gz
D3_chr21.vcf.gz.tbi	Vindija/VCF/Denisova/chr21_mq25_mapab100.vcf.gz.tbi
Chag8_chr21.vcf.gz	Chagyrskaya/VCF/chr21.noRB.vcf.gz
Chag8_chr21.vcf.gz.tbi	Chagyrskaya/VCF/chr21.noRB.vcf.gz.tbi
Mbuti_chr21.vcf.gz	altai/ModernHumans/vcf/SS6004471.hg19_1000g.21.mod.vcf.gz
Mbuti_chr21.vcf.gz.tbi	altai/ModernHumans/vcf/SS6004471.hg19_1000g.21.mod.vcf.gz.tbi
EOF
)

command -v curl >/dev/null || { echo "ERROR: hace falta curl" >&2; exit 1; }
command -v wget >/dev/null || { echo "ERROR: hace falta wget" >&2; exit 1; }

# --- elegir host ---
BASE=""
for h in "${HOSTS[@]}"; do
  if curl -sIf --max-time 20 "$h/Vindija/VCF/Altai/chr21_mq25_mapab100.vcf.gz" >/dev/null 2>&1; then
    BASE=$h; break
  fi
  echo "aviso: $h no responde, pruebo el siguiente" >&2
done
[ -n "$BASE" ] || { echo "ERROR: ningún servidor responde" >&2; exit 1; }
echo "servidor: $BASE"
echo

# --- comprobar las 10 URL antes de descargar ---
echo "Comprobando URL..."
missing=0
while IFS=$'\t' read -r local remote; do
  [ -n "$local" ] || continue
  size=$(curl -sIf --max-time 20 "$BASE/$remote" 2>/dev/null \
         | awk 'BEGIN{IGNORECASE=1}/^content-length:/{gsub(/\r/,"");print $2}' | tail -1)
  if [ -n "$size" ]; then
    printf "  ok     %-26s %10s bytes\n" "$local" "$size"
  else
    printf "  FALLA  %-26s  no encontrado: %s\n" "$local" "$remote"
    missing=$((missing+1))
  fi
done <<< "$FILES"

if [ "$missing" -gt 0 ]; then
  echo
  echo "ERROR: $missing fichero(s) no encontrados. No se descarga nada." >&2
  echo "Revisa las rutas en DATOS_NECESARIOS.md contra el índice del servidor." >&2
  exit 1
fi

# --- descargar ---
mkdir -p "$DEST"; cd "$DEST"
echo
echo "Descargando en $(pwd) ..."
while IFS=$'\t' read -r local remote; do
  [ -n "$local" ] || continue
  echo "--> $local"
  wget -c -q --show-progress -O "$local" "$BASE/$remote"
done <<< "$FILES"

# --- md5 de origen, para comprobar a mano ---
for d in Vindija/VCF/Altai Vindija/VCF/Vindija33.19 Vindija/VCF/Denisova Chagyrskaya/VCF; do
  curl -sf --max-time 30 -o "md5_$(echo "$d" | tr / _).txt" "$BASE/$d/md5sums.txt" || true
done

echo
echo "Descarga completa:"
ls -lh ./*.vcf.gz
echo
echo "Comprobación de integridad: los md5 de origen están en md5_*.txt. El nombre remoto es"
echo "chr21_mq25_mapab100.vcf.gz (o chr21.noRB.vcf.gz), el local lleva el del espécimen."
echo "Para comprobar uno:  md5sum D5_chr21.vcf.gz  y comparar con md5_Vindija_VCF_Altai.txt"
echo
echo "Comprobación rápida de que los VCF se abren:"
echo "  for f in *.vcf.gz; do echo -n \"\$f: \"; zcat \"\$f\" | head -1; done"
