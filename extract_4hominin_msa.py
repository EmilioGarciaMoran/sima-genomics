# RETIRADO (2026-10-21): NO USAR. Ver CORRECTION.md, sección "Fase 2".
# Este script usa accesiones escritas a mano que NO son de homínidos:
#   - KC758896 (Ust-Ishim) = Ostrinia scapulalis, transportador ABC (polilla)
# Cualquier resultado obtenido con él queda retirado.
from Bio import SeqIO

# Read alignment
alignment = SeqIO.parse("inputMitocondria/alignment_6_mtdna_CODING_aligned.fasta", "fasta")

# Keep only the 4 hominins (remove Pan NC_001643.1)
keep_ids = {
    "KF683087.1",   # Sima
    "FN673705.1",   # Denisova3
    "KC879692.1",   # Neandertal
    "KC758896.1",   # Ust-Ishim
    "NC_012920.1"   # Homo sapiens
}

records = [seq for seq in alignment if seq.id in keep_ids]

# Write clean 4-hominin alignment
SeqIO.write(records, "inputMitocondria/alignment_4hominin_CLEAN.fasta", "fasta")

print(f"✓ Extracted {len(records)} hominins (Pan removed)")
for r in records:
    print(f"  - {r.id} ({len(r.seq)} bp)")
