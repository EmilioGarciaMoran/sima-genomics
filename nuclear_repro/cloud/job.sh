#!/usr/bin/env bash
# Ejecuta UN trabajo (panel o vector de una ventana) con red de seguridad; lo llama run_vectors.sh.
#   - no empieza si ya pasó la hora límite (DEADLINE, epoch): queda pendiente para otra ejecución
#   - no empieza si necesita un panel que no existe (su ventana falló): queda como "sin_panel"
#   - tiempo máximo JOB_TIMEOUT segundos; si falla, se reintenta RETRIES veces
#   - siempre devuelve 0: un fallo no detiene los demás trabajos, solo se apunta en $STATE/failed
#
#   bash cloud/job.sh ETIQUETA PANEL_REQUERIDO|- 'comando'
set -uo pipefail
label=$1; needs=$2; cmd=$3
: "${STATE:?}" "${DEADLINE:?}" "${JOB_TIMEOUT:=3600}" "${RETRIES:=1}"

note() {  # note MENSAJE FICHERO — escritura con cerrojo: varios trabajos escriben a la vez
  { flock 9; printf '%s\t%s\t%s\n' "$(date -Is)" "$label" "$1" >> "$STATE/$2"; } 9> "$STATE/lock"
}

if [ "$(date +%s)" -ge "$DEADLINE" ]; then note "hora límite alcanzada" pending; exit 0; fi
if [ "$needs" != "-" ] && [ ! -s "$needs" ]; then note "sin panel (falló su ventana)" skipped; exit 0; fi

for try in $(seq 0 "$RETRIES"); do
  start=$(date +%s)
  out=$(timeout "$JOB_TIMEOUT" bash -c "$cmd" 2>&1)
  rc=$?
  if [ "$rc" = 0 ]; then
    note "$(( $(date +%s) - start )) s" done
    printf '%s\n' "$out" | tail -1
    exit 0
  fi
  reason=$([ "$rc" = 124 ] && echo "superó ${JOB_TIMEOUT} s" || echo "código $rc")
  echo "FALLO ($label, intento $((try + 1))): $reason" >&2
  printf '%s\n' "$out" | tail -3 | sed 's/^/    /' >&2
done
note "$reason: $(printf '%s' "$out" | tail -1 | tr '\t' ' ')" failed
exit 0
