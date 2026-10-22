#!/usr/bin/env bash
# NUBE — procesa las ventanas de pairwise/windows.bed: un panel y un vector por espécimen en cada una.
# Es acumulativo: lo que ya existe no se rehace. Añadir ventanas (cloud/make_windows.py) o especímenes
# (pairwise/specimens.tsv) y volver a lanzarlo solo calcula lo que falta.
# Al repositorio vuelve pairwise/windows/<id>/{panel.tsv.gz, vectors/*.vec.gz} (KB por ventana).
#
#   bash cloud/run_vectors.sh DATA_DIR REF.fa [hilos]
#   BUCKET=s3://sima-egarmo-2026 bash cloud/run_vectors.sh ...   # sube resultados y log a S3
#   PAIRWISE_DIR=otro/dir bash cloud/run_vectors.sh ...          # (los tests lo usan así)
#
# DATA_DIR contiene los BAM (indexados), VCF (bgzip + tabix) y FASTA de pairwise/specimens.tsv.
# Como solo se leen las ventanas, path también puede ser una URL (https://, s3://) a un BAM/VCF
# indexado: htslib descarga solo los tramos necesarios.
# Panel de cada ventana: sitios variables en los VCF de role=comparator y en el FASTA de role=outgroup.
# El panel de una ventana queda fijo una vez hecho: un comparador añadido después aporta su vector,
# pero no cambia qué sitios hay (bórralo a mano si quieres rehacerlo con el panel nuevo).
set -euo pipefail
cd "$(dirname "$0")/.."
DATA=${1:?DATA_DIR}; REF=${2:?REF.fa}; J=${3:-4}
PW=${PAIRWISE_DIR:-pairwise}
SPEC=$PW/specimens.tsv
WIN=$PW/windows.bed
[ -s "$WIN" ] || { echo "ERROR: falta $WIN (créalo con cloud/make_windows.py)" >&2; exit 1; }
mkdir -p "$PW/windows"
LOG=$PW/run_vectors.$(date +%Y%m%d_%H%M%S).log
exec > >(tee -a "$LOG") 2>&1
echo "inicio $(date -Is)  datos=$DATA  ref=$REF  hilos=$J"

rows() { grep -v '^#' "$SPEC" | awk -F'\t' 'NR>1 && NF>=6'; }
wins() { grep -v '^#' "$WIN" | awk -F'\t' 'NF>=4'; }
src() { case "$1" in *://*) echo "$1" ;; *) echo "$DATA/$1" ;; esac; }

vcfs=""; ogs=""
while IFS=$'\t' read -r name group role mode path sample; do
  [ "$role" = comparator ] && [ "$mode" = vcf ] && vcfs="$vcfs $(printf %q "$(src "$path")")"
  [ "$role" = outgroup ] && [ "$mode" = fasta ] && ogs="$ogs $(printf %q "$(src "$path")")"
done < <(rows)

# 1. Paneles que faltan (uno por ventana, en paralelo)
CMDS=$PW/windows/.panel_commands
wins | while IFS=$'\t' read -r chrom start end id; do
  [ -s "$PW/windows/$id/panel.tsv.gz" ] && continue
  printf 'python3 cloud/09_site_panel.py --ref %q --region %q --vcf %s --outgroup-fasta %s --out %q\n' \
    "$REF" "$chrom:$((start + 1))-$end" "$vcfs" "$ogs" "$PW/windows/$id/panel.tsv.gz"
done > "$CMDS"
echo "paneles por hacer: $(wc -l < "$CMDS")"
xargs -P "$J" -d '\n' -I{} bash -c {} < "$CMDS"

# 2. Vectores que faltan (ventana × espécimen, en paralelo)
CMDS=$PW/windows/.vector_commands
wins | while IFS=$'\t' read -r chrom start end id; do
  rows | while IFS=$'\t' read -r name group role mode path sample; do
    [ -s "$PW/windows/$id/vectors/$name.vec.gz" ] && continue
    extra=""; [ "$mode" = vcf ] && extra="--sample $(printf %q "$sample")"
    printf 'python3 cloud/10_specimen_vector.py --panel %q --name %q --%s %q %s --out %q\n' \
      "$PW/windows/$id/panel.tsv.gz" "$name" "$mode" "$(src "$path")" "$extra" "$PW/windows/$id/vectors"
  done
done > "$CMDS"
echo "vectores por hacer: $(wc -l < "$CMDS")"
xargs -P "$J" -d '\n' -I{} bash -c {} < "$CMDS"
rm -f "$PW"/windows/.*_commands

(cd "$PW" && find windows -name '*.gz' | sort | xargs md5sum > windows/MD5SUMS)
echo "fin $(date -Is)"
if [ -n "${BUCKET:-}" ]; then
  aws s3 sync "$PW/windows" "$BUCKET/pairwise/windows/"
  aws s3 cp "$WIN" "$BUCKET/pairwise/"
  aws s3 cp "$LOG" "$BUCKET/pairwise/logs/"
fi
