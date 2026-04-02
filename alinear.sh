#!/bin/bash

cd ~/sima/nuclear_analysis/bams

BWA="/usr/bin/bwa.kit/bwa"
SAMTOOLS="/usr/bin/samtools"
REF="$HOME/sima/ref/hg19.fa"
THREADS=8

echo "=== ALINEANDO A hg19 ==="

# DENISOVA 3 - ERR145625
$BWA mem -t $THREADS -M $REF \
  ~/sima/nuclear_analysis/fastqs/ERR145625/ERR145625_1.fastq.gz \
  ~/sima/nuclear_analysis/fastqs/ERR145625/ERR145625_2.fastq.gz 2>/dev/null | \
$SAMTOOLS view -b -S - | \
$SAMTOOLS sort -@ $THREADS -o denisova3_ERR145625.bam -
$SAMTOOLS index denisova3_ERR145625.bam
echo "✅ denisova3_ERR145625.bam"

# DENISOVA 3 - ERR145629
$BWA mem -t $THREADS -M $REF \
  ~/sima/nuclear_analysis/fastqs/ERR145629/ERR145629_1.fastq.gz \
  ~/sima/nuclear_analysis/fastqs/ERR145629/ERR145629_2.fastq.gz 2>/dev/null | \
$SAMTOOLS view -b -S - | \
$SAMTOOLS sort -@ $THREADS -o denisova3_ERR145629.bam -
$SAMTOOLS index denisova3_ERR145629.bam
echo "✅ denisova3_ERR145629.bam"

# DENISOVA 3 - ERR145630
$BWA mem -t $THREADS -M $REF \
  ~/sima/nuclear_analysis/fastqs/ERR145630/ERR145630_1.fastq.gz \
  ~/sima/nuclear_analysis/fastqs/ERR145630/ERR145630_2.fastq.gz 2>/dev/null | \
$SAMTOOLS view -b -S - | \
$SAMTOOLS sort -@ $THREADS -o denisova3_ERR145630.bam -
$SAMTOOLS index denisova3_ERR145630.bam
echo "✅ denisova3_ERR145630.bam"

# SIMA - ERR995360 (paired pero archivo único - usa 2 veces)
$BWA mem -t $THREADS -M $REF \
  ~/sima/nuclear_analysis/fastqs/ERR995360.fastq.gz \
  ~/sima/nuclear_analysis/fastqs/ERR995360.fastq.gz 2>/dev/null | \
$SAMTOOLS view -b -S - | \
$SAMTOOLS sort -@ $THREADS -o sima.bam -
$SAMTOOLS index sima.bam
echo "✅ sima.bam"

echo ""
echo "=== RESUMEN ==="
ls -lh *.bam
