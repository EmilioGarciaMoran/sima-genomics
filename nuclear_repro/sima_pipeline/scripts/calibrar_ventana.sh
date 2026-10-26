#!/usr/bin/env bash
# Calibración: ¿qué ventana de detección del daño usó Meyer 2016?
#
# El artículo dice "C->T terminal" sin precisar en cuántas posiciones. Con ventana 1
# obtenemos la mitad de las bases útiles que declara, y el cociente es muy constante entre
# los cinco especímenes (0,42 a 0,58), lo que apunta a un parámetro mal fijado y no a un
# error aleatorio.
#
# Esto NO es ajustar el análisis al resultado deseado: es determinar un parámetro que el
# artículo no especifica, usando como criterio el recuento de bases que él mismo publica.
# Queda documentado como tal, y las D-statistics se reportarán con todas las ventanas.
#
#   bash scripts/calibrar_ventana.sh [directorio_datos]
set -uo pipefail
cd "$(dirname "$0")/.."
DATA=${1:-$HOME/sima/sima_data}
REF="$DATA/ref/human_g1k_v37.fasta"
BAM="$DATA/bam"

echo "Calibración de la ventana de detección del daño"
echo "Criterio: acercarse a las bases útiles de Meyer 2016, Tabla 1"
echo
printf "%-20s %10s %12s %12s %9s\n" especimen ventana "bases útiles" Meyer cociente
while IFS=$'\t' read -r esp muestra corridas lecturas bases papel; do
  case "$esp" in \#*|especimen|"") continue ;; esac
  [ -s "$BAM/${esp}.bam" ] || continue
  for w in 1 2 3 4; do
    u=$(samtools view -h "$BAM/${esp}.bam" \
        | python3 scripts/damage_filter.py --ref "$REF" --library ss \
                  --damage-window "$w" --mask 3 2>&1 >/dev/null \
        | awk -F: '/BASES ÚTILES/{gsub(/[ ,]/,"",$2); print $2}')
    r=$(awk -v a="${u:-0}" -v b="$bases" 'BEGIN{if(b>0)printf "%.2f",a/b; else print "-"}')
    printf "%-20s %10s %12s %12s %9s\n" "$esp" "$w" "${u:-?}" "$bases" "$r"
  done
  echo
done < ESPECIMENES.tsv
echo "Se elige la ventana cuyo cociente quede más cerca de 1,00 de forma consistente"
echo "en los cinco especímenes."
