#!/usr/bin/env bash
# NUBE — construye el panel (si no existe) y un vector por espécimen de pairwise/specimens.tsv.
# Lo pesado (recorrer BAM y VCF completos) ocurre aquí; al repositorio solo vuelven
# pairwise/panel.tsv.gz y pairwise/vectors/*.vec.gz (pocos MB).
#
#   bash cloud/run_vectors.sh DATA_DIR REF.fa [hilos]
#   BUCKET=s3://sima-egarmo-2026 bash cloud/run_vectors.sh ...   # sube resultados y log a S3
#   PAIRWISE_DIR=otro/dir CHROMS="chr1 chr2" bash cloud/run_vectors.sh ...   # (los tests lo usan así)
#
# DATA_DIR contiene los BAM (indexados), VCF (bgzip + tabix) y FASTA de pairwise/specimens.tsv.
# Paneles: sitios variables en los VCF de role=comparator y en el FASTA de role=outgroup.
set -euo pipefail
cd "$(dirname "$0")/.."
DATA=${1:?DATA_DIR}; REF=${2:?REF.fa}; J=${3:-4}
PW=${PAIRWISE_DIR:-pairwise}
CHROMS=${CHROMS:-$(printf 'chr%s ' {1..22})}
SPEC=$PW/specimens.tsv
LOG=$PW/run_vectors.$(date +%Y%m%d_%H%M%S).log
mkdir -p "$PW/vectors"
exec > >(tee -a "$LOG") 2>&1
echo "inicio $(date -Is)  datos=$DATA  ref=$REF  hilos=$J"

rows() { grep -v '^#' "$SPEC" | awk -F'\t' 'NR>1 && NF>=6'; }

if [ ! -s "$PW/panel.tsv.gz" ]; then
  vcfs=$(rows | awk -F'\t' -v d="$DATA" '$3=="comparator" && $4=="vcf"{printf "%s/%s ", d, $5}')
  ogs=$(rows | awk -F'\t' -v d="$DATA" '$3=="outgroup" && $4=="fasta"{printf "%s/%s ", d, $5}')
  python3 cloud/09_site_panel.py --ref "$REF" --vcf $vcfs --outgroup-fasta $ogs \
    --chroms $CHROMS --out "$PW/panel.tsv.gz"
else
  echo "$PW/panel.tsv.gz ya existe: se reutiliza (bórralo para regenerarlo)"
fi

# Un proceso por espécimen, J a la vez; si alguno falla, xargs devuelve error y el script se detiene
CMDS=$PW/vectors/.commands
rows | while IFS=$'\t' read -r name group role mode path sample; do
  extra=""; [ "$mode" = vcf ] && extra="--sample $(printf %q "$sample")"
  printf 'python3 cloud/10_specimen_vector.py --panel %q --name %q --%s %q %s --out %q\n' \
    "$PW/panel.tsv.gz" "$name" "$mode" "$DATA/$path" "$extra" "$PW/vectors"
done > "$CMDS"
xargs -P "$J" -d '\n' -I{} bash -c {} < "$CMDS"
rm -f "$CMDS"

(cd "$PW" && md5sum panel.tsv.gz vectors/*.vec.gz > vectors/MD5SUMS)
echo "fin $(date -Is)"
if [ -n "${BUCKET:-}" ]; then
  aws s3 cp "$PW/panel.tsv.gz" "$BUCKET/pairwise/"
  aws s3 sync "$PW/vectors" "$BUCKET/pairwise/vectors/"
  aws s3 cp "$LOG" "$BUCKET/pairwise/logs/"
fi
