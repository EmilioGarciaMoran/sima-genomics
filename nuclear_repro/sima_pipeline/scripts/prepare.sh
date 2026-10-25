#!/usr/bin/env bash
# Preparación del análisis de Sima: referencia, índice de bwa, metadatos de ENA y FASTQ.
# Todo es trabajo preparatorio, sin ninguna decisión analítica: se puede dejar corriendo.
#
#   bash scripts/prepare.sh [directorio_datos]
#   nohup bash scripts/prepare.sh > results/preparacion.log 2>&1 &
#
# Unas 2 h, dominadas por `bwa index` (~90 min a un hilo, no paralelizable).
# Disco: ~9 GB (referencia 3 GB + índice 5 GB + FASTQ ~0,2 GB).
#
# REANUDABLE: salta cada paso cuyo resultado ya exista.
set -uo pipefail
cd "$(dirname "$0")/.."

DATA=${1:-$HOME/sima/sima_data}
REF="$DATA/ref/human_g1k_v37.fasta"
FQ="$DATA/fastq"
META="$DATA/meta"
mkdir -p "$DATA/ref" "$FQ" "$META" results

echo "=================================================================="
echo "Preparación del análisis de Sima  |  inicio: $(date '+%F %T')"
echo "  datos en: $DATA"
echo "=================================================================="

# --- 0. herramientas -----------------------------------------------------
falta=0
for t in bwa samtools wget curl; do
  if command -v "$t" >/dev/null; then
    printf "  ok     %-10s %s\n" "$t" "$(command -v $t)"
  else
    printf "  FALTA  %-10s\n" "$t"; falta=1
  fi
done
if [ "$falta" -ne 0 ]; then
  echo
  echo "Instala lo que falte y vuelve a lanzarlo:" >&2
  echo "  sudo apt install bwa samtools wget curl" >&2
  exit 1
fi

L=$(df -BG --output=avail "$DATA" 2>/dev/null | tail -1 | tr -dc '0-9')
echo "  disco libre: ${L} GB (hacen falta ~9)"
[ -n "$L" ] && [ "$L" -lt 12 ] && { echo "ERROR: poco espacio" >&2; exit 1; }

# --- 1. metadatos de ENA: qué espécimen es cada corrida -------------------
# Pendiente desde octubre y decisivo: todo el diseño exige analizar cada espécimen
# por separado, y aún no sabemos cuál es cuál.
echo
echo "--- 1. Metadatos de ENA (PRJEB10597) ---"
MF="$META/PRJEB10597.tsv"
if [ -s "$MF" ]; then
  echo "  ya descargado"
else
  curl -sf --max-time 60 -o "$MF" \
    "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB10597&result=read_run&fields=run_accession,sample_accession,sample_title,sample_alias,library_name,description,read_count,base_count,fastq_ftp,fastq_bytes&format=tsv" \
    || { echo "  AVISO: no se pudo descargar el filereport (se reintenta en la próxima ejecución)" >&2; rm -f "$MF"; }
fi
if [ -s "$MF" ]; then
  echo "  corridas encontradas: $(($(wc -l < "$MF") - 1))"
  echo "  --- contenido (columnas clave) ---"
  awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}
              {printf "   %-14s %-14s %-28s %12s lecturas\n",
                      $h["run_accession"], $h["sample_accession"],
                      substr($h["sample_title"]" "$h["sample_alias"],1,28), $h["read_count"]}' "$MF"
fi

# --- 2. referencia GRCh37 -------------------------------------------------
# Se usa human_g1k_v37, cuyos cromosomas se llaman "21" y no "chr21", igual que los
# VCF de MPI-EVA. Mezclar convenciones de nombres arruinaría el cruce de posiciones.
echo
echo "--- 2. Referencia GRCh37 ---"
# Por HTTPS: el FTP del EBI va lento y wget -c sobre él puede dejar bytes de más.
URL=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
if [ -s "$REF" ]; then
  echo "  ya descomprimida: $(du -h "$REF" | cut -f1)"
else
  # Sin -c: reanudar sobre un fichero previo puede concatenar en vez de continuar,
  # y deja "trailing garbage" al final del flujo gzip.
  rm -f "$REF.gz"
  echo "  descargando (~850 MB)..."
  wget -q -O "$REF.gz" "$URL" || {
    echo "  ERROR: no se pudo descargar la referencia" >&2
    echo "  Pruébalo a mano: wget -O $REF.gz $URL" >&2; exit 1; }
  echo "  descomprimiendo..."
  # gunzip avisa de bytes sobrantes tras el flujo y devuelve error, pero el contenido
  # descomprimido es correcto. La verificación de verdad son las longitudes de abajo.
  gunzip -c "$REF.gz" > "$REF" 2>/dev/null
  if [ ! -s "$REF" ]; then
    echo "  ERROR: la descompresión no produjo nada" >&2; rm -f "$REF"; exit 1
  fi
  rm -f "$REF.gz"
fi
[ -s "$REF.fai" ] || { echo "  indexando con samtools faidx..."; samtools faidx "$REF" || {
  echo "  ERROR: samtools faidx falló; la referencia está incompleta" >&2
  rm -f "$REF" "$REF.fai"; exit 1; }; }

# VERIFICACIÓN REAL: longitudes de los 22 autosomas en GRCh37. Si cuadran, la referencia
# está completa y es la correcta. Es más fiable que cualquier prueba sobre el .gz.
echo "  verificando longitudes de los autosomas..."
python3 - "$REF.fai" <<'PYEOF' || { echo "  ERROR: la referencia no es GRCh37 o está incompleta" >&2; exit 1; }
import sys
ESPERADO = {"1":249250621,"2":243199373,"3":198022430,"4":191154276,"5":180915260,
            "6":171115067,"7":159138663,"8":146364022,"9":141213431,"10":135534747,
            "11":135006516,"12":133851895,"13":115169878,"14":107349540,"15":102531392,
            "16":90354753,"17":81195210,"18":78077248,"19":59128983,"20":63025520,
            "21":48129895,"22":51304566}
visto = {}
for linea in open(sys.argv[1]):
    c = linea.split("\t")
    if c[0] in ESPERADO:
        visto[c[0]] = int(c[1])
malos = [c for c in ESPERADO if visto.get(c) != ESPERADO[c]]
if malos:
    print("    DISCREPAN:", ", ".join(sorted(malos, key=int)))
    for c in sorted(malos, key=int)[:5]:
        print(f"      cr{c}: {visto.get(c, 'ausente')} en vez de {ESPERADO[c]}")
    sys.exit(1)
print(f"    los 22 autosomas cuadran con GRCh37")
PYEOF
echo "  cromosomas: $(cut -f1 "$REF.fai" | head -5 | tr '\n' ' ')..."
echo "  (sin prefijo chr, como los VCF de MPI-EVA)"

# --- 3. índice de bwa -----------------------------------------------------
echo
echo "--- 3. Índice de bwa (el paso largo, ~90 min) ---"
if [ -s "$REF.bwt" ] && [ -s "$REF.sa" ]; then
  echo "  ya existe"
else
  echo "  empezando a las $(date '+%T')..."
  bwa index "$REF" 2>&1 | sed 's/^/    /' || { echo "  ERROR en bwa index" >&2; exit 1; }
  echo "  terminado a las $(date '+%T')"
fi

# --- 4. FASTQ de Sima -----------------------------------------------------
echo
echo "--- 4. FASTQ de Sima desde ENA ---"
if [ -s "$MF" ]; then
  n=0
  while IFS=$'\t' read -r run ftp; do
    [ -n "$run" ] && [ "$run" != "run_accession" ] || continue
    IFS=';' read -ra urls <<< "$ftp"
    for u in "${urls[@]}"; do
      [ -n "$u" ] || continue
      f="$FQ/$(basename "$u")"
      if [ -s "$f" ]; then echo "  ya está: $(basename "$f")"; continue; fi
      wget -c -q -O "$f" "https://$u" && gzip -t "$f" 2>/dev/null \
        && echo "  ok: $(basename "$f") ($(du -h "$f" | cut -f1))" \
        || { echo "  FALLO: $u" >&2; rm -f "$f"; }
    done
    n=$((n+1))
  done < <(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}
                       {print $h["run_accession"]"\t"$h["fastq_ftp"]}' "$MF")
  echo "  corridas procesadas: $n"
  echo "  total descargado: $(du -sh "$FQ" | cut -f1)"
else
  echo "  SALTADO: hacen falta los metadatos de ENA (paso 1)"
fi

echo
echo "=================================================================="
echo "Preparación terminada: $(date '+%F %T')"
echo "  referencia: $REF"
echo "  índice:     $([ -s "$REF.bwt" ] && echo 'listo' || echo 'FALTA')"
echo "  FASTQ:      $(ls "$FQ" 2>/dev/null | wc -l) ficheros en $FQ"
echo "  metadatos:  $MF"
echo
echo "Siguiente paso (mañana): mapeo con bwa aln sin siembra, filtro de daño"
echo "y extracción de alelos. Ese código se escribe con sus pruebas, no deprisa."
echo "=================================================================="
