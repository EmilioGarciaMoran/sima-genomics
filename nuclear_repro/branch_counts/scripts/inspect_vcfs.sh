#!/usr/bin/env bash
# Inspecciona la cabecera y los primeros registros de cada VCF descargado, para saber
# qué campos hay antes de escribir el lector. No carga los ficheros enteros: se corta
# en cuanto tiene lo que necesita.
#
#   bash scripts/inspect_vcfs.sh ~/sima/public_vcfs/chr21
#
# Pega la salida completa: es corta.
set -uo pipefail
DIR=${1:-data/chr21}
cd "$DIR" || { echo "ERROR: no existe $DIR" >&2; exit 1; }

echo "=== Ficheros ==="
ls -lh ./*.vcf.gz
echo
echo "=== Comprobación md5 frente al origen ==="
chk() {  # chk <fichero_local> <fichero_md5> <nombre_remoto>
  if [ -f "$2" ]; then
    want=$(awk -v n="$3" '$2==n || $2=="./"n {print $1}' "$2" | head -1)
    if [ -n "$want" ]; then
      got=$(md5sum "$1" | cut -d' ' -f1)
      [ "$got" = "$want" ] && echo "  ok     $1" || echo "  DIFIERE $1  (local $got, origen $want)"
    else
      echo "  ?      $1: no encuentro $3 dentro de $2"
    fi
  else
    echo "  ?      $1: falta $2"
  fi
}
chk D5_chr21.vcf.gz     md5_Vindija_VCF_Altai.txt        chr21_mq25_mapab100.vcf.gz
chk Vi3319_chr21.vcf.gz md5_Vindija_VCF_Vindija33.19.txt chr21_mq25_mapab100.vcf.gz
chk D3_chr21.vcf.gz     md5_Vindija_VCF_Denisova.txt     chr21_mq25_mapab100.vcf.gz
chk Chag8_chr21.vcf.gz  md5_Chagyrskaya_VCF.txt          chr21.noRB.vcf.gz
echo "  (Mbuti no trae md5sums.txt en su directorio)"

for f in D5 Vi3319 D3 Chag8 Mbuti; do
  v="${f}_chr21.vcf.gz"
  [ -f "$v" ] || { echo; echo "=== $v: NO ESTÁ ==="; continue; }
  echo
  echo "================================================================"
  echo "=== $v"
  echo "================================================================"
  echo "--- versión, referencia y contig:"
  zcat "$v" 2>/dev/null | head -400 | grep -E '^##(fileformat|reference|source|contig=<ID=(21|chr21)[,>])' | head -6
  echo "--- campos INFO declarados:"
  zcat "$v" 2>/dev/null | head -400 | grep '^##INFO' | sed 's/^##INFO=<//; s/>$//' | cut -c1-110
  echo "--- campos FORMAT declarados:"
  zcat "$v" 2>/dev/null | head -400 | grep '^##FORMAT' | sed 's/^##FORMAT=<//; s/>$//' | cut -c1-110
  echo "--- FILTER declarados:"
  zcat "$v" 2>/dev/null | head -400 | grep '^##FILTER' | sed 's/^##FILTER=<//; s/>$//' | cut -c1-110
  echo "--- línea de columnas (nombres de muestra):"
  zcat "$v" 2>/dev/null | grep -m1 '^#CHROM'
  echo "--- 3 primeros registros (recortados a 220 caracteres):"
  zcat "$v" 2>/dev/null | grep -v '^#' | head -3 | cut -c1-220
  echo "--- ¿hay alelo ancestral declarado? (AA, AAncestral, ANC...):"
  zcat "$v" 2>/dev/null | head -400 | grep -iE '^##INFO.*ID=(AA|AAncestral|ANC|ANCESTRAL)' || echo "  (no aparece)"
done

echo
echo "================================================================"
echo "=== Densidad de registros: primeras 200.000 líneas de datos de D5"
echo "================================================================"
zcat D5_chr21.vcf.gz 2>/dev/null | grep -v '^#' | head -200000 | \
  awk '{n++; if(n==1)a=$2; b=$2} END{printf "  %d registros entre la posicion %d y %d (%.0f kb)\n", n, a, b, (b-a)/1000}'
echo
echo "=== Reparto de FILTER en esas 200.000 líneas (D5):"
zcat D5_chr21.vcf.gz 2>/dev/null | grep -v '^#' | head -200000 | cut -f7 | sort | uniq -c | sort -rn | head -10
echo
echo "Listo. Pega toda la salida."
