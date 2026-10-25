#!/usr/bin/env bash
# Segunda inspección: registros de variante reales y campos ancestrales del VCF de Mbuti.
# Necesario para escribir el lector sin suponer nada.
#
#   bash scripts/inspect_vcfs2.sh ~/sima/public_vcfs/chr21
#
# Lee solo el principio de cada fichero. Pega la salida completa.
set -uo pipefail
DIR=${1:-data/chr21}
cd "$DIR" || { echo "ERROR: no existe $DIR" >&2; exit 1; }
N=${2:-3000000}     # líneas de datos a recorrer por fichero

echo "=== md5 de Chagyrskaya: buscando el listado en el servidor ==="
for u in https://ftp.eva.mpg.de/neandertal/Chagyrskaya/VCF/md5sums.txt \
         https://ftp.eva.mpg.de/neandertal/Chagyrskaya/VCF/md5sum.txt \
         https://ftp.eva.mpg.de/neandertal/Chagyrskaya/md5sums.txt; do
  if curl -sIf --max-time 15 "$u" >/dev/null 2>&1; then
    echo "  existe: $u"; curl -sf --max-time 20 "$u" | grep -i chr21 || true
  else
    echo "  no existe: $u"
  fi
done
echo "  md5 local de Chag8_chr21.vcf.gz:"
md5sum Chag8_chr21.vcf.gz

for f in D5 Vi3319 D3 Chag8; do
  v="${f}_chr21.vcf.gz"
  echo
  echo "================================================================"
  echo "=== $v : registros con ALT distinto de '.' (sitios variables)"
  echo "================================================================"
  zcat "$v" 2>/dev/null | grep -v '^#' | head -"$N" | awk -F'\t' '$5!="."' | head -6 | cut -c1-200
  echo "--- reparto de genotipo (campo GT) en $N líneas:"
  zcat "$v" 2>/dev/null | grep -v '^#' | head -"$N" | cut -f10 | cut -d: -f1 | sort | uniq -c | sort -rn | head -8
  echo "--- reparto de ALT en $N líneas:"
  zcat "$v" 2>/dev/null | grep -v '^#' | head -"$N" | cut -f5 | sort | uniq -c | sort -rn | head -8
done

echo
echo "================================================================"
echo "=== Mbuti: campos ancestrales"
echo "================================================================"
echo "--- 4 registros CON CAnc declarado:"
zcat Mbuti_chr21.vcf.gz 2>/dev/null | grep -v '^#' | head -"$N" | grep -m4 'CAnc=' | cut -c1-260
echo
echo "--- 3 registros con TSseq (bases de los 6 primates):"
zcat Mbuti_chr21.vcf.gz 2>/dev/null | grep -v '^#' | head -"$N" | grep -m3 'TSseq=' | cut -f1,2,4,5,8 | cut -c1-260
echo
echo "--- valores que toma TS (orden de especies):"
zcat Mbuti_chr21.vcf.gz 2>/dev/null | grep -v '^#' | head -"$N" | grep -o 'TS=[^;]*' | sort | uniq -c | sort -rn | head -5
echo
echo "--- 4 sitios variables en Mbuti (ALT distinto de '.'):"
zcat Mbuti_chr21.vcf.gz 2>/dev/null | grep -v '^#' | head -"$N" | awk -F'\t' '$5!="."' | head -4 | cut -c1-260
echo
echo "--- de $N lineas, cuantas traen cada campo:"
zcat Mbuti_chr21.vcf.gz 2>/dev/null | grep -v '^#' | head -"$N" | \
  awk -F'\t' '{n++; if($8~/CAnc=/)c++; if($8~/TSseq=/)t++; if($8~/Map20=/)m++; if($8~/(^|;)RM(;|$)/)r++}
              END{printf "  total %d | CAnc %d | TSseq %d | Map20 %d | RM %d\n", n,c,t,m,r}'

echo
echo "================================================================"
echo "=== Solape de posiciones entre D5 y Mbuti (primeras $N lineas de cada uno)"
echo "================================================================"
zcat D5_chr21.vcf.gz    2>/dev/null | grep -v '^#' | head -"$N" | cut -f2 | sort -n > /tmp/_p_d5.txt
zcat Mbuti_chr21.vcf.gz 2>/dev/null | grep -v '^#' | head -"$N" | cut -f2 | sort -n > /tmp/_p_mb.txt
echo -n "  D5:    "; wc -l < /tmp/_p_d5.txt
echo -n "  Mbuti: "; wc -l < /tmp/_p_mb.txt
echo -n "  comunes: "; comm -12 /tmp/_p_d5.txt /tmp/_p_mb.txt | wc -l
echo "  rango D5:    $(head -1 /tmp/_p_d5.txt) - $(tail -1 /tmp/_p_d5.txt)"
echo "  rango Mbuti: $(head -1 /tmp/_p_mb.txt) - $(tail -1 /tmp/_p_mb.txt)"
rm -f /tmp/_p_d5.txt /tmp/_p_mb.txt
echo
echo "Listo. Pega toda la salida."
