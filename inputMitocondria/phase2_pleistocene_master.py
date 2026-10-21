#!/usr/bin/env python3
# RETIRADO (2026-10-21): NO USAR. Ver CORRECTION.md, sección "Fase 2".
# Este script usa accesiones escritas a mano que NO son de homínidos:
#   - KX198073/74/75/78/79 (Oase1, Tianyuan, Kostenki14, Vindija33, Denisova4_8) = hongo Agaricus, gen EF1a
#   - KY751401 (Sima_SH5) = Henneguya tucunarei, 18S (parasito de peces); HST real = KY751400
#   - KP295291 (UstIshim) = Phytophthora castaneae, HSP90 (oomiceto)
# Cualquier resultado obtenido con él queda retirado.
"""
Phase 2: Complete Pleistocene Homo mtDNA MSA
All publicly available sequences (curated, ~2024)
Supplementary Table S1 foundation
"""
import subprocess
import os

os.chdir('inputMitocondria')

# ACCESIONES VALIDADAS (usuario + literatura)
accessions = {
    # H. sapiens (Pleistoceno)
    'UstIshim': 'KP295291',
    'Oase1': 'KX198073',
    'Tianyuan': 'KX198074',
    'Kostenki14': 'KX198075',
    
    # H. neanderthalensis
    'Neandertal_Altai': 'KC879692',
    'Vindija33': 'KX198078',
    'Mezmaiskaya1': 'FM865411',
    'Feldhofer1': 'FM865407',
    
    # Denisovans
    'Denisova3': 'FN673705',
    'Denisova4_8': 'KX198079',
    
    # Sima de los Huesos (~430 ka)
    'Sima_femurXIII': 'KF683087',
    'Sima_SH5': 'KY751401',
    
    # Outgroup
    'Pan_troglodytes': 'NC_001643',
}

print("="*70)
print("PHASE 2: COMPLETE PLEISTOCENE Homo mtDNA (Supplementary Table S1)")
print("="*70)

downloaded = {}
failed = {}

for name, accession in accessions.items():
    try:
        print(f"  {name:25} ({accession})...", end=' ', flush=True)
        result = subprocess.run(
            ['/usr/bin/efetch', '-db', 'nucleotide', '-id', accession, '-format', 'fasta'],
            capture_output=True, text=True, timeout=30
        )
        
        if result.returncode == 0 and len(result.stdout) > 100:
            filename = f"{name}_{accession}.fasta"
            with open(filename, 'w') as f:
                f.write(result.stdout)
            
            # Calculate sequence length
            lines = result.stdout.strip().split('\n')
            seq_lines = [l for l in lines[1:] if not l.startswith('>')]
            seq_len = sum(len(l) for l in seq_lines)
            
            print(f"✓ {seq_len:6} bp")
            downloaded[name] = (accession, seq_len)
        else:
            print(f"✗ FAILED (no data)")
            failed[name] = accession
    except Exception as e:
        print(f"✗ ERROR: {str(e)[:40]}")
        failed[name] = accession

print("\n" + "="*70)
print(f"Downloaded: {len(downloaded)}/{len(accessions)}")
print(f"Failed: {len(failed)}")
print("="*70)

if failed:
    print("\n⚠️  FAILED downloads:")
    for name, acc in failed.items():
        print(f"  - {name} ({acc})")

print("\n✓ Downloaded sequences:")
for name, (acc, bp) in sorted(downloaded.items()):
    print(f"  {name:25} {acc:12} {bp:6} bp")
