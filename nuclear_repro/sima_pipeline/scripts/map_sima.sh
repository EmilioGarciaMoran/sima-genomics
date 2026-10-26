#!/usr/bin/env bash
# Pieza 1 del pipeline: mapeo de cada espécimen de Sima.
#
#   bash scripts/map_sima.sh [directorio_datos] [hilos]
#   nohup bash scripts/map_sima.sh > results/mapeo.log 2>&1 &
#
# Parámetros fijados en ../PREREGISTRO_D_statistics.md, apartado 10 ter. No se cambian
# sin documentar la desviación con su motivo y su fecha.
#
# Un BAM por espécimen, nunca fusionados entre sí. Las corridas de una MISMA muestra sí se
# juntan: son librerías del mismo hueso (ver ESPECIMENES.tsv).
# La escápula AT-6672 se mapea igualmente, para poder enseñar su contaminación, pero queda
# excluida del análisis.
#
# Unos 15 min en total. Reanudable: salta el espécimen cuyo BAM ya exista.
set -uo pipefail
cd "$(dirname "$0")/.."

DATA=${1:-$HOME/sima/sima_data}
T=${2:-4}
REF="$DATA/ref/human_g1k_v37.fasta"
FQ="$DATA/fastq"
OUT="$DATA/bam"
TMP="$DATA/_tmp_map"
mkdir -p "$OUT" "$TMP" results

# --- parámetros preregistrados -------------------------------------------
ALN_L=16500      # siembra desactivada: la siembra penaliza el daño del extremo 5'
ALN_N=0.01       # desajustes permitidos
ALN_O=2          # huecos permitidos
MIN_MAPQ=30      # literal de Meyer 2016: "a map quality score of 30 or greater"
MIN_LEN=35       # por debajo, el mapeo único no es fiable

echo "=================================================================="
echo "Mapeo de Sima  |  inicio: $(date '+%F %T')"
echo "  referencia: $REF"
echo "  parámetros: bwa aln -l $ALN_L -n $ALN_N -o $ALN_O | mapq>=$MIN_MAPQ | len>=$MIN_LEN"
echo "  salida: $OUT"
echo "=================================================================="

for t in bwa samtools; do
  command -v "$t" >/dev/null || { echo "ERROR: falta $t (sudo apt install bwa samtools)" >&2; exit 1; }
done
[ -s "$REF.bwt" ] || { echo "ERROR: falta el índice de bwa. Ejecuta scripts/prepare.sh" >&2; exit 1; }

# --- un espécimen por fila de ESPECIMENES.tsv ----------------------------
total_ok=0; fallos=""
while IFS=$'\t' read -r esp muestra corridas lecturas bases papel; do
  case "$esp" in \#*|especimen|"") continue ;; esac

  BAM="$OUT/${esp}.bam"
  if [ -s "$BAM" ]; then
    echo "[$esp] ya mapeado, se salta"
    total_ok=$((total_ok+1)); continue
  fi

  echo
  echo "------------------------------------------------------------------"
  echo "[$esp] $(date '+%T')  muestra $muestra, papel: $papel"

  # juntar las corridas de ESTA muestra (y solo de esta)
  CAT="$TMP/${esp}.fastq.gz"
  rm -f "$CAT"
  faltan=""
  for run in ${corridas//,/ }; do
    f="$FQ/${run}.fastq.gz"
    [ -s "$f" ] || { faltan="$faltan $run"; continue; }
    cat "$f" >> "$CAT"
  done
  if [ -n "$faltan" ]; then
    echo "[$esp] FALLO: faltan los FASTQ de:$faltan" >&2
    fallos="$fallos $esp"; rm -f "$CAT"; continue
  fi
  echo "[$esp] corridas juntadas: $corridas  ($(du -h "$CAT" | cut -f1))"

  # mapeo, filtrado y eliminación de duplicados en una sola tubería
  SAI="$TMP/${esp}.sai"
  if ! bwa aln -t "$T" -l "$ALN_L" -n "$ALN_N" -o "$ALN_O" "$REF" "$CAT" > "$SAI" 2> "$TMP/${esp}.aln.log"; then
    echo "[$esp] FALLO en bwa aln (ver $TMP/${esp}.aln.log)" >&2
    fallos="$fallos $esp"; continue
  fi
  if ! bwa samse "$REF" "$SAI" "$CAT" 2>> "$TMP/${esp}.aln.log" \
       | python3 scripts/filter_dedup.py --min-mapq "$MIN_MAPQ" --min-len "$MIN_LEN" \
                 --etiqueta "$esp" \
       | samtools sort -@ "$T" -T "$TMP/sort_$esp" -o "$BAM" - 2>/dev/null; then
    echo "[$esp] FALLO en samse/filtrado/ordenado" >&2
    rm -f "$BAM"; fallos="$fallos $esp"; continue
  fi
  samtools index "$BAM"

  n=$(samtools view -c "$BAM")
  echo "[$esp] BAM: $(du -h "$BAM" | cut -f1), $n lecturas conservadas"
  echo "[$esp] lecturas de partida según ENA: $lecturas  -> conservadas: $n"
  rm -f "$SAI" "$CAT"
  total_ok=$((total_ok+1))
done < ESPECIMENES.tsv

echo
echo "=================================================================="
echo "Fin: $(date '+%F %T')   especímenes mapeados: $total_ok"
[ -n "$fallos" ] && echo "FALLARON:$fallos  (relanza el script: es reanudable)"
echo
echo "Resumen:"
printf "  %-22s %12s %12s %8s\n" especimen conservadas ENA pct
for b in "$OUT"/*.bam; do
  [ -s "$b" ] || continue
  e=$(basename "$b" .bam)
  n=$(samtools view -c "$b")
  ena=$(awk -F'\t' -v e="$e" '$1==e{print $4}' ESPECIMENES.tsv)
  pct=$(awk -v a="$n" -v b="${ena:-0}" 'BEGIN{if(b>0)printf "%.1f%%", 100*a/b; else print "-"}')
  printf "  %-22s %12s %12s %8s\n" "$e" "$n" "${ena:--}" "$pct"
done
echo
echo "Siguiente pieza: filtro de daño (C->T terminal) y enmascarado 3+3."
echo "=================================================================="
