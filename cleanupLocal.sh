#!/bin/bash

# =============================================================================
# LIMPIEZA LOCAL - Borrar FASTQ parciales + verificar espacio
# =============================================================================

echo "🧹 LIMPIEZA LOCAL DE ARCHIVOS PARCIALES"
echo "========================================"
echo ""

# Directorio de trabajo
WORKDIR="$HOME/sima"

# Espacio ANTES
echo "📊 ESPACIO EN DISCO ANTES:"
df -h "$WORKDIR" | tail -1
echo ""

# Ver qué vamos a borrar
echo "📂 Archivos a borrar:"
echo ""

# FASTQ parciales (si existen)
if [ -d "$WORKDIR/nuclear_analysis/fastqs" ]; then
    echo "En nuclear_analysis/fastqs/:"
    du -sh "$WORKDIR/nuclear_analysis/fastqs"/* 2>/dev/null || echo "  (vacío)"
    echo ""
fi

if [ -d "$WORKDIR/aws/fastqs" ]; then
    echo "En aws/fastqs/:"
    du -sh "$WORKDIR/aws/fastqs"/* 2>/dev/null || echo "  (vacío)"
    echo ""
fi

# Calcular espacio total a liberar
SPACE_TO_FREE=$(du -sh "$WORKDIR/nuclear_analysis/fastqs" "$WORKDIR/aws/fastqs" 2>/dev/null | awk '{sum+=$1} END {print sum}')

echo "💾 Espacio a liberar: ~$SPACE_TO_FREE"
echo ""

read -p "¿Borrar archivos FASTQ parciales locales? (s/n): " CONFIRM

if [ "$CONFIRM" != "s" ]; then
    echo "Cancelado"
    exit 0
fi

# BORRAR
echo ""
echo "🗑️  Borrando archivos..."

# Borrar FASTQ parciales
rm -rf "$WORKDIR/nuclear_analysis/fastqs/ERR"*.fastq.gz 2>/dev/null
rm -rf "$WORKDIR/nuclear_analysis/fastqs/SRR"*.fastq.gz 2>/dev/null
rm -rf "$WORKDIR/aws/fastqs" 2>/dev/null

# Borrar BAMs vacíos/incorrectos
rm -f "$WORKDIR/nuclear_analysis/bams/Sima_aligned.bam" 2>/dev/null
rm -f "$WORKDIR/nuclear_analysis/bams/Ust-Ishim_aligned.bam" 2>/dev/null
rm -f "$WORKDIR/nuclear_analysis/bams/Neandertal_aligned.bam" 2>/dev/null
rm -f "$WORKDIR/nuclear_analysis/bams/"*.bai 2>/dev/null

# Borrar logs temporales
rm -f "$WORKDIR/nuclear_analysis/pipeline_"*.log 2>/dev/null
rm -f "$WORKDIR/nuclear_analysis/"*.pid 2>/dev/null

echo "✅ Archivos borrados"
echo ""

# Espacio DESPUÉS
echo "📊 ESPACIO EN DISCO DESPUÉS:"
df -h "$WORKDIR" | tail -1
echo ""

# Detalle por partición
echo "╔════════════════════════════════════════════════════════════╗"
echo "║           ESPACIO EN DISCO - TODAS LAS PARTICIONES          ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""
df -h | grep -E "Filesystem|/$|/home|/mnt"
echo ""

# Lo que queda en nuclear_analysis
echo "📂 CONTENIDO RESTANTE EN ~/sima/nuclear_analysis/:"
echo ""
du -sh "$WORKDIR/nuclear_analysis"/* 2>/dev/null | sort -h || echo "  (vacío)"
echo ""

# Resumen
echo "╔════════════════════════════════════════════════════════════╗"
echo "║                    LIMPIEZA COMPLETA                        ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""
echo "✅ Archivos parciales borrados"
echo "✅ BAMs incorrectos borrados"
echo "✅ Logs temporales borrados"
echo ""
echo "🎯 Sistema listo para pipeline AWS (descarga directa en nube)"
echo ""
