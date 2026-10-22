#!/usr/bin/env bash
# Paso 3b (opcional): las mismas verosimilitudes con ANGSD y los parámetros exactos del manuscrito
#   -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -minMapQ 20 -minQ 20, sitios = CDS fusionadas,
# más una segunda ejecución solo con transversiones (-noTrans 1).
# Se salta sin error si angsd no está instalado (p. ej. en la integración continua).
#
#   bash scripts/03b_angsd.sh [slices] [catalog] [results]
set -euo pipefail
cd "$(dirname "$0")/.."
SLICES=${1:-slices}; CATALOG=${2:-catalog}; OUT=${3:-results}/angsd

if ! command -v angsd >/dev/null; then
  echo "  angsd no está instalado: se omite el paso 3b"
  exit 0
fi
mkdir -p "$OUT"
angsd -version 2>&1 | head -1 > "$OUT/VERSION.txt" || true

# Fichero de sitios (1-based) y de regiones a partir del BED de CDS
grep -v '^#' "$CATALOG/cds_regions.bed" | awk 'NF>=3{for(p=$2+1;p<=$3;p++) print $1"\t"p}' \
  | sort -k1,1 -k2,2n -u > "$OUT/cds.sites"
cut -f1 "$OUT/cds.sites" | uniq > "$OUT/cds.chrs"
angsd sites index "$OUT/cds.sites" 2>/dev/null

run() {  # $1 = etiqueta, $2 = lista de BAM
  for nt in 0 1; do
    suf=$([ "$nt" = 1 ] && echo ".notrans" || echo "")
    angsd -bam "$2" -sites "$OUT/cds.sites" -rf "$OUT/cds.chrs" \
      -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -noTrans "$nt" \
      -minMapQ 20 -minQ 20 -out "$OUT/$1$suf" 2> "$OUT/$1$suf.log"
    n=$( (zcat "$OUT/$1$suf.mafs.gz" 2>/dev/null || true) | tail -n +2 | wc -l)
    echo "  ANGSD $1$suf: $n sitios"
  done
}

ls "$SLICES"/*.bam > "$OUT/POOL.bamlist"
while read -r b; do
  s=$(basename "$b" .bam); echo "$b" > "$OUT/$s.bamlist"; run "$s" "$OUT/$s.bamlist"
done < "$OUT/POOL.bamlist"
run POOL "$OUT/POOL.bamlist"
