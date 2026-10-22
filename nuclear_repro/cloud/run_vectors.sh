#!/usr/bin/env bash
# NUBE — procesa las ventanas pendientes de pairwise/windows.bed: un panel y un vector por espécimen.
# Pensado para lanzarlo a mano y que termine solo sin sorpresas de coste ni pérdidas:
#
#   PLAN=1      bash cloud/run_vectors.sh DATA_DIR REF.fa       # 1) ver qué haría, sin calcular nada
#   MAX_WINDOWS=1 BUCKET=s3://... bash cloud/run_vectors.sh ...  # 2) piloto: una sola ventana
#   BUCKET=s3://... bash cloud/run_vectors.sh DATA_DIR REF.fa 8 # 3) el resto
#
# Variables (todas opcionales):
#   BUCKET          s3://... destino; se comprueba que se puede escribir ANTES de calcular nada
#   MAX_WINDOWS     procesar como mucho N ventanas pendientes en esta ejecución (0 = todas)
#   MAX_HOURS       no empezar trabajos nuevos pasadas N horas (por defecto 4); lo no empezado
#                   queda pendiente para la siguiente ejecución
#   JOB_TIMEOUT     segundos máximos por trabajo (por defecto 3600); RETRIES reintentos (1)
#   HEARTBEAT_MIN   minutos entre latidos de progreso en S3 (por defecto 5)
#   SHUTDOWN=1      apagar la instancia al terminar (también si algo falla)
#   PAIRWISE_DIR    directorio de trabajo (por defecto pairwise; los tests usan otro)
#
# Pase lo que pase (fallos, hora límite, Ctrl-C), al salir: resumen, MD5SUMS y subida a S3 de todo lo
# terminado. Un trabajo que falla no detiene a los demás; queda en run_*/failed y se reintenta en la
# siguiente ejecución, porque solo se da por hecho lo que existe. Código de salida 1 si hubo fallos.
#
# DATA_DIR contiene los BAM (indexados), VCF (bgzip + tabix) y FASTA de pairwise/specimens.tsv;
# path también puede ser una URL a un fichero indexado (htslib lee solo los tramos de la ventana).
# Panel de cada ventana: sitios variables en los VCF de role=comparator y el FASTA de role=outgroup;
# queda fijo una vez hecho (un comparador añadido después solo aporta su vector).
set -euo pipefail
cd "$(dirname "$0")/.."
DATA=${1:?DATA_DIR}; REF=${2:?REF.fa}; J=${3:-4}
PW=${PAIRWISE_DIR:-pairwise}
SPEC=$PW/specimens.tsv; WIN=$PW/windows.bed
MAX_WINDOWS=${MAX_WINDOWS:-0}; MAX_HOURS=${MAX_HOURS:-4}
export JOB_TIMEOUT=${JOB_TIMEOUT:-3600} RETRIES=${RETRIES:-1}
HEARTBEAT_MIN=${HEARTBEAT_MIN:-5}; SHUTDOWN=${SHUTDOWN:-0}; PLAN=${PLAN:-0}; BUCKET=${BUCKET:-}

die() { echo "ERROR: $*" >&2; exit 1; }
[ -s "$WIN" ] || die "falta $WIN (créalo con cloud/make_windows.py)"
[ -s "$SPEC" ] || die "falta $SPEC"
[ -s "$REF" ] || die "no existe la referencia $REF"
for c in python3 timeout flock; do command -v $c >/dev/null || die "falta $c"; done
python3 -c "import pysam, numpy" 2>/dev/null || die "python3 sin pysam/numpy (micromamba activate nucrepro)"

rows() { grep -v '^#' "$SPEC" | awk -F'\t' 'NR>1 && NF>=6'; }
wins() { grep -v '^#' "$WIN" | awk -F'\t' 'NF>=4'; }
src() { case "$1" in *://*) echo "$1" ;; *) echo "$DATA/$1" ;; esac; }

# Una sola ejecución a la vez sobre el mismo directorio (dos a la vez repetirían trabajo y coste)
mkdir -p "$PW"
exec 8> "$PW/.run_lock"
flock -n 8 || die "ya hay otra ejecución en marcha sobre $PW (cerrojo $PW/.run_lock)"

# Temporales de una ejecución cortada a medias (nunca cuentan como hechos; se rehacen)
find "$PW" -name '.*.tmp' -delete 2>/dev/null || true

# Comprobar que todos los ficheros de entrada existen antes de empezar (las URL no se comprueban)
missing=""
while IFS=$'\t' read -r name group role mode path sample; do
  case "$path" in *://*) continue ;; esac
  [ -e "$DATA/$path" ] || missing="$missing\n  $name: $DATA/$path"
done < <(rows)
[ -z "$missing" ] || die "faltan ficheros de entrada:$(printf "$missing")"

# Ventanas pendientes (sin panel o con algún vector por hacer), en el orden de windows.bed
mapfile -t names < <(rows | cut -f1)
pending=()
while IFS=$'\t' read -r chrom start end id; do
  todo=0; [ -s "$PW/windows/$id/panel.tsv.gz" ] || todo=1
  for n in "${names[@]}"; do [ -s "$PW/windows/$id/vectors/$n.vec.gz" ] || todo=1; done
  [ "$todo" = 1 ] && pending+=("$chrom"$'\t'"$start"$'\t'"$end"$'\t'"$id")
done < <(wins)
total_pending=${#pending[@]}
[ "$MAX_WINDOWS" -gt 0 ] && [ "$total_pending" -gt "$MAX_WINDOWS" ] && pending=("${pending[@]:0:$MAX_WINDOWS}")

vcfs=""; ogs=""
while IFS=$'\t' read -r name group role mode path sample; do
  [ "$role" = comparator ] && [ "$mode" = vcf ] && vcfs="$vcfs $(printf %q "$(src "$path")")"
  [ "$role" = outgroup ] && [ "$mode" = fasta ] && ogs="$ogs $(printf %q "$(src "$path")")"
done < <(rows)

TS=$(date +%Y%m%d_%H%M%S)
STATE=$PW/runs/$TS
PANELS=$(mktemp); VECTORS=$(mktemp)
for w in "${pending[@]}"; do
  IFS=$'\t' read -r chrom start end id <<< "$w"
  panel=$PW/windows/$id/panel.tsv.gz
  if [ ! -s "$panel" ]; then
    printf 'bash cloud/job.sh %q - %q\n' "$id/panel" \
      "python3 cloud/09_site_panel.py --ref $(printf %q "$REF") --region $chrom:$((start + 1))-$end --vcf $vcfs --outgroup-fasta $ogs --out $(printf %q "$panel")" >> "$PANELS"
  fi
  while IFS=$'\t' read -r name group role mode path sample; do
    [ -s "$PW/windows/$id/vectors/$name.vec.gz" ] && continue
    extra=""; [ "$mode" = vcf ] && extra="--sample $(printf %q "$sample")"
    printf 'bash cloud/job.sh %q %q %q\n' "$id/$name" "$panel" \
      "python3 cloud/10_specimen_vector.py --panel $(printf %q "$panel") --name $(printf %q "$name") --$mode $(printf %q "$(src "$path")") $extra --out $(printf %q "$PW/windows/$id/vectors")" >> "$VECTORS"
  done < <(rows)
done
n_panels=$(wc -l < "$PANELS"); n_vectors=$(wc -l < "$VECTORS")

echo "== Plan"
echo "  ventanas pendientes: $total_pending; en esta ejecución: ${#pending[@]}$([ "$MAX_WINDOWS" -gt 0 ] && echo " (MAX_WINDOWS=$MAX_WINDOWS)")"
echo "  trabajos: $n_panels paneles + $n_vectors vectores, $J a la vez"
echo "  límites: $JOB_TIMEOUT s por trabajo, $RETRIES reintento(s), no se empieza nada pasadas $MAX_HOURS h"
echo "  S3: ${BUCKET:-(sin BUCKET: los resultados quedan solo en esta máquina)}; apagado al terminar: $([ "$SHUTDOWN" = 1 ] && echo sí || echo no)"
if [ "$PLAN" = 1 ]; then rm -f "$PANELS" "$VECTORS"; echo "(PLAN=1: no se ejecuta nada)"; exit 0; fi
if [ $((n_panels + n_vectors)) = 0 ]; then rm -f "$PANELS" "$VECTORS"; echo "Nada pendiente."; exit 0; fi

# Comprobar S3 antes de gastar una hora de cálculo
if [ -n "$BUCKET" ]; then
  command -v aws >/dev/null || die "BUCKET definido pero falta aws"
  echo "prueba $(date -Is)" > /tmp/.s3_probe_$$
  aws s3 cp /tmp/.s3_probe_$$ "$BUCKET/pairwise/runs/.probe" >/dev/null || die "no se puede escribir en $BUCKET"
  rm -f /tmp/.s3_probe_$$
fi

mkdir -p "$STATE"; touch "$STATE"/{done,failed,pending,skipped}
cp "$PANELS" "$STATE/panel_jobs"; cp "$VECTORS" "$STATE/vector_jobs"; rm -f "$PANELS" "$VECTORS"
export STATE DEADLINE=$(( $(date +%s) + MAX_HOURS * 3600 ))
# el tee ignora Ctrl-C: así sigue vivo para registrar el resumen y la subida tras una interrupción
exec > >(trap '' INT TERM; exec tee -a "$STATE/log.txt") 2>&1
echo "inicio $(date -Is)  datos=$DATA  ref=$REF  estado=$STATE"

progress() { echo "$(wc -l < "$STATE/done") hechos, $(wc -l < "$STATE/failed") fallidos, de $((n_panels + n_vectors))"; }

heartbeat() {
  while sleep $((HEARTBEAT_MIN * 60)); do
    echo "$(date -Is)  $(progress)" > "$STATE/HEARTBEAT.txt"
    [ -n "$BUCKET" ] && aws s3 cp "$STATE/HEARTBEAT.txt" "$BUCKET/pairwise/runs/$TS/" >/dev/null 2>&1
  done
}

finish() {
  rc=$?
  set +e   # aquí nada puede abortar: lo importante es que la subida a S3 ocurra siempre
  trap - EXIT
  trap '' INT TERM PIPE   # un segundo Ctrl-C tampoco interrumpe la subida
  if [ -n "${HB:-}" ]; then pkill -P "$HB" 2>/dev/null; kill "$HB" 2>/dev/null; fi
  # lo que no llegó a lanzarse (Ctrl-C, error del propio script) también queda como pendiente
  (cd "$PW" && find windows -name '*.gz' ! -name '.*' | sort | xargs -r md5sum > windows/MD5SUMS) || true
  nf=$(wc -l < "$STATE/failed")
  {
    echo "Ejecución $TS: $(date -Is)"
    echo "  $(progress)"
    echo "  sin empezar por la hora límite: $(wc -l < "$STATE/pending"); sin panel: $(wc -l < "$STATE/skipped")"
    echo "  ventanas aún pendientes en windows.bed: $(( total_pending - ${#pending[@]} )) + las de esta ejecución que no se completaron"
    if [ "$nf" -gt 0 ]; then echo "  FALLIDOS (se reintentan en la siguiente ejecución):"; cut -f2,3 "$STATE/failed" | sed 's/^/    /'; fi
  } | tee "$STATE/SUMMARY.txt"
  if [ -n "$BUCKET" ]; then
    echo "subiendo a $BUCKET ..."
    for t in 1 2 3; do
      aws s3 sync "$PW/windows" "$BUCKET/pairwise/windows/" --exclude '.*' --exclude '*/.*' \
        && aws s3 sync "$STATE" "$BUCKET/pairwise/runs/$TS/" && aws s3 cp "$WIN" "$BUCKET/pairwise/" && break
      echo "  la subida falló (intento $t); reintento en 30 s"; sleep 30
    done || true
  fi
  if [ "$SHUTDOWN" = 1 ]; then echo "apagando la instancia"; sudo shutdown -h now; fi
  if [ "$nf" -gt 0 ]; then exit 1; fi
  exit "$rc"
}
trap finish EXIT
trap 'echo "interrumpido"; exit 130' INT TERM

heartbeat > /dev/null 2>&1 & HB=$!   # sin heredar la salida: su sleep no retiene el registro
echo "== Paneles ($n_panels)"
xargs -P "$J" -d '\n' -I{} bash -c {} < "$STATE/panel_jobs"
echo "== Vectores ($n_vectors)"
xargs -P "$J" -d '\n' -I{} bash -c {} < "$STATE/vector_jobs"
echo "fin del cálculo $(date -Is)"
