#!/usr/bin/env bash
# Paso 5: diagnóstico de los resultados ORIGINALES de la fase 3a (marzo 2026), sin re-ejecutar nada.
# Aplica el mismo recuento de biparticiones del paso 3 a los ficheros que hay en
# inputMitocondria/phase3a_results/, para documentar qué decían realmente.
#
#   bash scripts/05_legacy_check.sh
#
# Salida: results/LEGACY_PHASE3A_CHECK.txt
# Ejecutar desde mtDNA_repro/.
set -euo pipefail
L=../inputMitocondria/phase3a_results
O=results/LEGACY_PHASE3A_CHECK.txt
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT

# Los ficheros de RAxML legacy usan nombres largos; se recortan a las 10 letras de los de IQ-TREE
long2short='s/Denisova3_70ka/Denisova3_/g; s/Neandertal_Vindija_44ka/Neandertal/g;
            s/Feldhofer1_40ka/Feldhofer1/g; s/Mezmaiskaya1_60ka/Mezmaiskay/g; s/Pan_outgroup/Pan_outgro/g'
sed "$long2short" "$L/raxml_gtr.raxml.bootstraps" > "$TMP/rx_boot.nwk"

{
  echo "== ¿Se enraizó el árbol IQ-TREE original? =="
  grep -m1 "UNROOTED" "$L/iqtree_hky.iqtree" || echo "(sin aviso)"
  echo
  echo "== Árbol de consenso original (iqtree_hky.contree) =="
  cat "$L/iqtree_hky.contree"
  echo
  echo "== Agrupaciones presentes en el árbol ML original (100 = presente, 0 = ausente) =="
  python3 -I scripts/clade_support.py --header
  python3 -I scripts/clade_support.py "$L/iqtree_hky.treefile" "phase3a_IQTREE_HKY+G4_ML_tree"
  echo
  echo "== Réplicas bootstrap originales de RAxML-NG (raxml_gtr.raxml.bootstraps) =="
  python3 -I scripts/clade_support.py --header
  python3 -I scripts/clade_support.py "$TMP/rx_boot.nwk" "phase3a_RAxML_GTR+G_stdboot"
  echo
  echo "Lectura: el nodo con soporte 93 del árbol original separa {Sima, Denisova} de"
  echo "{Neandertales, Pan}. Enraizado en Pan, ese nodo es el clado Sima+Denisova."
} | tee "$O"
