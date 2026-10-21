#!/usr/bin/env bash
# Paso 1: obtener las 6 secuencias mitocondriales sin alinear.
#
#   bash scripts/01_get_sequences.sh genbank   # descarga de NCBI (recomendado, en tu máquina)
#   bash scripts/01_get_sequences.sh recover   # quita los huecos del alineamiento legacy de phase3a
#
# En modo genbank se conserva el título original de GenBank en la cabecera (title=...).
# Al final se ejecuta scripts/validate_sequences.py, que aborta si algo no cuadra.
# Salida: sequences/<taxon>.fasta, sequences/all_raw.fasta, sequences/SEQUENCES.tsv (longitud, nº de N, sha256)
# Ejecutar desde mtDNA_repro/.
set -euo pipefail
MODE=${1:?usa: genbank | recover}
RAW=sequences
LEGACY=../inputMitocondria/phase3a_results/alignment_6magnificos.phy
mkdir -p "$RAW"

case "$MODE" in
  genbank)
    tail -n +2 accessions.tsv | while IFS=$'\t' read -r taxon acc _; do
      curl -sSf "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${acc}&rettype=fasta&retmode=text" \
        | sed "1s/^>[^ ]* \(.*\)$/>${taxon} ${acc} source=GenBank title=\1/" > "$RAW/${taxon}.fasta"
      sleep 0.4
    done ;;
  recover)
    python3 -I - "$LEGACY" "$RAW" accessions.tsv <<'PY'
import sys
phy, raw, acc_tsv = sys.argv[1:]
acc = {l.split('\t')[0]: l.split('\t')[1] for l in open(acc_tsv).read().splitlines()[1:]}
for line in open(phy).read().splitlines()[1:]:
    if not line.strip(): continue
    name, seq = line[:10].strip(), line[10:]
    ungapped = seq.replace('-', '').upper()
    with open(f"{raw}/{name}.fasta", "w") as f:
        f.write(f">{name} {acc[name]} source=degapped_from_phase3a_alignment\n")
        for i in range(0, len(ungapped), 70): f.write(ungapped[i:i+70] + "\n")
PY
    ;;
  *) echo "modo desconocido: $MODE" >&2; exit 1 ;;
esac

# Unir en el orden del manifiesto y documentar
: > "$RAW/all_raw.fasta"
printf "taxon\taccession\tlength\tN_count\tsha256_sequence\n" > "$RAW/SEQUENCES.tsv"
tail -n +2 accessions.tsv | while IFS=$'\t' read -r taxon acc _; do
  cat "$RAW/${taxon}.fasta" >> "$RAW/all_raw.fasta"
  seq=$(grep -v '>' "$RAW/${taxon}.fasta" | tr -d '\n' | tr 'acgtn' 'ACGTN')
  n=$(printf "%s" "$seq" | tr -cd 'N' | wc -c)
  sha=$(printf "%s" "$seq" | sha256sum | cut -c1-16)
  printf "%s\t%s\t%s\t%s\t%s\n" "$taxon" "$acc" "${#seq}" "$n" "$sha" >> "$RAW/SEQUENCES.tsv"
done
cat "$RAW/SEQUENCES.tsv"

# Validación obligatoria antes de alinear (aborta si alguna secuencia no es lo que dice ser)
python3 -I scripts/validate_sequences.py accessions.tsv "$RAW"
