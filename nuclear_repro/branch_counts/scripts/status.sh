#!/usr/bin/env bash
# Estado del recorrido genómico en una pantalla.
#   bash scripts/status.sh            una vez
#   watch -n 30 bash scripts/status.sh   refrescando
cd "$(dirname "$0")/.."
WORK=${1:-$HOME/sima/public_vcfs/_tmp}
N=$(ls results/genome/chr*_inf.npz 2>/dev/null | wc -l)
echo "== $(date '+%F %T') =="
if pgrep -f run_genome.sh >/dev/null; then
  echo "estado:      EN MARCHA (pid $(pgrep -f run_genome.sh | head -1))"
else
  echo "estado:      PARADO"
fi
echo "cromosomas:  $N de 22 terminados"
if [ "$(pgrep -cf run_genome.sh)" -gt 1 ] 2>/dev/null; then
  echo "AVISO:       hay MÁS DE UN proceso en marcha. Pára todos con: pkill -f run_genome.sh"
fi
if [ "$N" -gt 0 ]; then
  echo -n "             hechos: "
  ls results/genome/chr*_inf.npz 2>/dev/null | sed 's/.*chr//; s/_inf.npz//' | sort -n | tr '\n' ' '
  echo
  echo "             sitios acumulados: $(python3 - <<'PY' 2>/dev/null || echo "n/d"
import glob, numpy as np
print(f"{sum(len(np.load(f)['pos']) for f in glob.glob('results/genome/chr*_inf.npz')):,}")
PY
)"
fi
echo "en curso:    $(du -sh "$WORK" 2>/dev/null | cut -f1 || echo '-') descargados en $WORK"
echo "disco libre: $(df -h --output=avail "$WORK" 2>/dev/null | tail -1 | tr -d ' ')"
echo "salidas:     $(du -sh results/genome 2>/dev/null | cut -f1)"
echo "-- últimas líneas del registro --"
tail -5 results/genoma.log 2>/dev/null
