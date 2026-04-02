#!/usr/bin/env python3
"""
Phase 2: Download 6 Magnificent mtDNA
CORRECTED ACCESSIONS (validated: DeepSeek + ChatGPT)
"""
import subprocess
import os

os.chdir('inputMitocondria')

# VALIDATED ACCESSIONS (mtDNA complete genomes ONLY)
accessions = {
    'Sima_femurXIII_430ka': 'KF683087',        # ✓ CORRECT (Meyer et al., 2014)
    'Denisova3_70ka': 'NC_013993',             # ✓ CORRECTED (was FN673705 - nuclear)
    'Neandertal_Vindija_44ka': 'NC_020627',   # ✓ CORRECTED (was KC879692 - SRA)
    'Feldhofer1_40ka': 'FM865407',             # ✓ CORRECT (Green et al., 2008)
    'Mezmaiskaya1_60ka': 'FM865411',           # ✓ CORRECT (Green et al., 2008)
    'Pan_troglodytes': 'NC_001643',           # ✓ CORRECT (reference)
}

print("="*80)
print("PHASE 2: DOWNLOAD 6 MAGNIFICENT mtDNA (CORRECTED ACCESSIONS)")
print("="*80)
print("\nValidation source: DeepSeek + ChatGPT editorial review")
print("Standard: mtDNA complete genomes ONLY (~16.5 kb each)\n")

downloaded = {}
failed = {}

for name, accession in accessions.items():
    try:
        print(f"  {name:35} ({accession:12})...", end=' ', flush=True)
        result = subprocess.run(
            ['/usr/bin/efetch', '-db', 'nucleotide', '-id', accession, '-format', 'fasta'],
            capture_output=True, text=True, timeout=30
        )
        
        if result.returncode == 0 and len(result.stdout) > 100:
            filename = f"{name}_{accession}.fasta"
            with open(filename, 'w') as f:
                f.write(result.stdout)
            
            lines = result.stdout.strip().split('\n')
            seq_lines = [l for l in lines[1:] if not l.startswith('>')]
            seq_len = sum(len(l) for l in seq_lines)
            
            print(f"✓ {seq_len:6} bp")
            downloaded[name] = (accession, seq_len)
        else:
            print(f"✗ FAILED")
            failed[name] = accession
    except Exception as e:
        print(f"✗ ERROR: {str(e)[:30]}")
        failed[name] = accession

print("\n" + "="*80)
print(f"Downloaded: {len(downloaded)}/{len(accessions)}")
print("="*80)

if downloaded:
    print("\n✓ SEQUENCES DOWNLOADED:")
    total_bp = 0
    for name, (acc, bp) in sorted(downloaded.items()):
        status = "COMPLETE" if bp >= 16500 else "PARTIAL"
        print(f"  {name:35} {acc:12} {bp:6} bp  [{status}]")
        total_bp += bp
    print(f"\n  Total bp: {total_bp}")

if failed:
    print("\n⚠️  FAILED DOWNLOADS:")
    for name, acc in failed.items():
        print(f"  {name} ({acc})")

print("\n" + "="*80)
print("Next: run phase2_compile_CORRECTED.py")
print("="*80)
