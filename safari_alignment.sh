#!/bin/bash

set -e

WORK_DIR="$HOME/sima"
SAFARI_DIR="$WORK_DIR/SAFARI"
ANALYSIS_DIR="$WORK_DIR/nuclear_analysis"
FASTQ_DIR="$ANALYSIS_DIR/fastqs"
BAM_DIR="$ANALYSIS_DIR/bams_safari"
REF="$WORK_DIR/ref/hg19.fa"

echo "════════════════════════════════════════════════════════════"
echo "ALINEAMIENTO - Denisova 3 + Sima madre"
echo "════════════════════════════════════════════════════════════"
echo ""

mkdir -p $BAM_DIR
cd $BAM_DIR

echo "[1] Samples a alinear:"
RUNS=(
  "ERR145625:Denisova3_run1"
  "ERR145629:Denisova3_run2"
  "ERR145630:Denisova3_run3"
  "ERR995360:Sima_madre"
)

for run_info in "${RUNS[@]}"; do
  RUN=$(echo $run_info | cut -d: -f1)
  LABEL=$(echo $run_info | cut -d: -f2)
  if [ -f "$FASTQ_DIR/$RUN/${RUN}_1.fastq.gz" ]; then
    SIZE=$(du -sh "$FASTQ_DIR/$RUN/" | cut -f1)
    echo "   ✅ $LABEL: $SIZE"
  fi
done
echo ""

align_bwa() {
  local RUN=$1
  local LABEL=$2
  
  echo "───────────────────────────────────────────────────────"
  echo "ALINEANDO: $LABEL ($RUN)"
  echo "───────────────────────────────────────────────────────"
  
  local R1="$FASTQ_DIR/$RUN/${RUN}_1.fastq.gz"
  local R2="$FASTQ_DIR/$RUN/${RUN}_2.fastq.gz"
  
  time (
    bwa aln -n 0.01 -o 2 -l 16500 -t 4 $REF $R1 > ${LABEL}_1.sai
    bwa aln -n 0.01 -o 2 -l 16500 -t 4 $REF $R2 > ${LABEL}_2.sai
    bwa sampe -P $REF ${LABEL}_1.sai ${LABEL}_2.sai $R1 $R2 | \
      samtools view -bS - | \
      samtools sort -o ${LABEL}.bam -
  ) 2>&1 | tee ${LABEL}.log
  
  samtools index ${LABEL}.bam
  rm -f ${LABEL}_*.sai
  
  echo "Estadísticas:"
  samtools flagstat ${LABEL}.bam
  echo ""
}

echo "[2] Iniciando alineamientos..."
START=$(date +%s)

for run_info in "${RUNS[@]}"; do
  RUN=$(echo $run_info | cut -d: -f1)
  LABEL=$(echo $run_info | cut -d: -f2)
  
  if [ -f "$FASTQ_DIR/$RUN/${RUN}_1.fastq.gz" ]; then
    align_bwa "$RUN" "$LABEL"
  fi
done

END=$(date +%s)
ELAPSED=$((END - START))

echo "════════════════════════════════════════════════════════════"
echo "✅ ALINEAMIENTOS COMPLETADOS"
echo "════════════════════════════════════════════════════════════"
echo ""
echo "BAMs generados:"
ls -lh $BAM_DIR/*.bam
echo ""
echo "Tiempo: $(printf '%dh %dm %ds' $((ELAPSED/3600)) $((ELAPSED%3600/60)) $((ELAPSED%60)))"
echo ""

ls -1 $BAM_DIR/*.bam > $BAM_DIR/bam_list.txt
echo "✅ BAM list guardado"
