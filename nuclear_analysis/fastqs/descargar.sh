#!/bin/bash
# DESCARGAS NUCLEARES - HOMINIDOS

cd ~/sima/nuclear_analysis/fastqs/

echo "=== DENISOVA 3 (ERP001519) ==="
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR145/623/ERR145623/ERR145623_1.fastq.gz &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR145/623/ERR145623/ERR145623_2.fastq.gz &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR145/625/ERR145625/ERR145625_1.fastq.gz &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR145/625/ERR145625/ERR145625_2.fastq.gz &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR145/627/ERR145627/ERR145627_1.fastq.gz &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR145/627/ERR145627/ERR145627_2.fastq.gz &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR145/628/ERR145628/ERR145628_1.fastq.gz &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR145/628/ERR145628/ERR145628_2.fastq.gz &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR145/629/ERR145629/ERR145629_1.fastq.gz &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR145/629/ERR145629/ERR145629_2.fastq.gz &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR145/630/ERR145630/ERR145630_1.fastq.gz &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR145/630/ERR145630/ERR145630_2.fastq.gz &

echo "=== SIMA (PRJEB10597) ==="
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR995/360/ERR995360/ERR995360_1.fastq.gz &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR995/360/ERR995360/ERR995360_2.fastq.gz &

echo "=== UST-ISHIM (PRJEB6622) ==="
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR566/093/ERR566093/ERR566093_1.fastq.gz &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR566/093/ERR566093/ERR566093_2.fastq.gz &

echo "=== NEANDERTAL (PRJEB12635) ==="
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR229/913/ERR229913/ERR229913_1.fastq.gz &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR229/913/ERR229913/ERR229913_2.fastq.gz &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR229/910/ERR229910/ERR229910_1.fastq.gz &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR229/910/ERR229910/ERR229910_2.fastq.gz &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR229/912/ERR229912/ERR229912_1.fastq.gz &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR229/912/ERR229912/ERR229912_2.fastq.gz &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR229/911/ERR229911/ERR229911_1.fastq.gz &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR229/911/ERR229911/ERR229911_2.fastq.gz &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR229/914/ERR229914/ERR229914_1.fastq.gz &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR229/914/ERR229914/ERR229914_2.fastq.gz &

wait
echo "=== ✅ DESCARGAS COMPLETAS ==="
