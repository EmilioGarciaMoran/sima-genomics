#!/bin/bash

set -e

WORK_DIR="$HOME/sima"
SAFARI_DIR="$WORK_DIR/SAFARI"

echo "════════════════════════════════════════════════════════════"
echo "SAFARI SETUP - Iniciando..."
echo "════════════════════════════════════════════════════════════"
echo ""

# [1] CREAR ESTRUCTURA
mkdir -p $SAFARI_DIR/{graph,binaries}
cd $SAFARI_DIR

# [2] DESCARGAR PANGENOMA (DTU)
echo "[1/4] Descargando pangenoma hominin de DTU..."
cd $SAFARI_DIR/graph
wget -nc -r -l1 --no-parent -nH --cut-dirs=2 \
  ftp://ftp.healthtech.dtu.dk:/public/SAFARI_graph/
echo "✅ Pangenoma descargado"

# [3] DESCARGAR BINARIO
echo ""
echo "[2/4] Descargando binario SAFARI..."
cd $SAFARI_DIR/binaries
wget -c "https://github.com/grenaud/SAFARI/releases/download/zenodo_release/safari" \
  -O safari 2>/dev/null || echo "Intento Zenodo..."
[ -f safari ] && chmod +x safari && echo "✅ Binario descargado"

# [4] VERIFICACIÓN
echo ""
echo "[3/4] Verificación..."
ls -lh $SAFARI_DIR/graph/ | grep -E "\.og|\.gbwt|\.min|\.ry" | head -5
ls -lh $SAFARI_DIR/binaries/

echo ""
echo "════════════════════════════════════════════════════════════"
echo "✅ SAFARI SETUP COMPLETADO"
echo "════════════════════════════════════════════════════════════"
echo ""
echo "FASTQs listos para alineamiento:"
ls -lh ~/sima/nuclear_analysis/fastqs/ERR*/ERR*_1.fastq.gz
echo ""
echo "Pangenoma + índices RYmer descargados en:"
echo "  $SAFARI_DIR/graph/"
echo ""
echo "Binario SAFARI en:"
echo "  $SAFARI_DIR/binaries/safari"
echo ""
