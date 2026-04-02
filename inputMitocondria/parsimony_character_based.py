#!/usr/bin/env python3
"""
CHARACTER-BASED PARSIMONY VALIDATION
Test H1 vs H2 topologies using actual nucleotide characters
Ancestral state determined by Pan outgroup

Usage: python3 parsimony_character_based.py
"""

from Bio import SeqIO
from collections import defaultdict

print("="*80)
print("CHARACTER-BASED PARSIMONY ANALYSIS (with Ancestral State Determination)")
print("="*80)

# Read alignment
alignment_file = 'alignment_6magnificos.fasta'
sequences = {}
for record in SeqIO.parse(alignment_file, 'fasta'):
    sequences[record.id] = str(record.seq).upper()

print(f"\nSequences loaded: {len(sequences)}")
for name, seq in sequences.items():
    print(f"  {name:30s} {len(seq):6d} bp")

# Get alignment length
align_length = len(list(sequences.values())[0])
print(f"\nAlignment length: {align_length} bp")

# Define taxa
taxa_order = [
    'Pan_outgroup',
    'Sima_430ka',
    'Denisova3_70ka',
    'Neandertal_Vindija_44ka',
    'Feldhofer1_40ka',
    'Mezmaiskaya1_60ka',
]

# Verify all taxa present
missing = [t for t in taxa_order if t not in sequences]
if missing:
    print(f"ERROR: Missing taxa: {missing}")
    exit(1)

# ============================================================================
# STEP 1: Extract parsimony-informative sites
# ============================================================================
print("\n" + "="*80)
print("STEP 1: IDENTIFY PARSIMONY-INFORMATIVE SITES")
print("="*80)

pars_info_sites = []
variable_sites = []
parsimony_informative_chars = []

for pos in range(align_length):
    # Get characters at this position
    chars = {}
    for taxon in taxa_order:
        c = sequences[taxon][pos]
        if c not in ['-', 'N']:  # Exclude gaps and ambiguities
            chars[taxon] = c
    
    # Skip if position is not fully covered
    if len(chars) < len(taxa_order):
        continue
    
    # Check if position is variable
    unique_states = set(chars.values())
    if len(unique_states) == 1:
        continue  # Not variable
    
    variable_sites.append((pos, chars))
    
    # Check if parsimony-informative (≥2 taxa with state A, ≥2 taxa with state B)
    state_counts = defaultdict(int)
    for state in unique_states:
        count = sum(1 for c in chars.values() if c == state)
        state_counts[state] = count
    
    # Parsimony-informative if at least 2 states with ≥2 taxa each
    informative_states = sum(1 for count in state_counts.values() if count >= 2)
    if informative_states >= 2:
        pars_info_sites.append((pos, chars))
        parsimony_informative_chars.append(chars)

print(f"\nTotal positions: {align_length}")
print(f"Variable sites: {len(variable_sites)}")
print(f"Parsimony-informative sites: {len(pars_info_sites)}")
print(f"Percentage informative: {100*len(pars_info_sites)/align_length:.2f}%")

# ============================================================================
# STEP 2: Determine ancestral state for each site
# ============================================================================
print("\n" + "="*80)
print("STEP 2: ANCESTRAL STATE DETERMINATION (Pan as outgroup)")
print("="*80)

ancestral_states = {}

for pos, chars in pars_info_sites:
    pan_state = chars['Pan_outgroup']
    ancestral_states[pos] = pan_state

print(f"\nAncestral state determined for {len(ancestral_states)} sites")
print("(Pan state = ancestral; alternative states = derived)")

# ============================================================================
# STEP 3: Classify derived state changes
# ============================================================================
print("\n" + "="*80)
print("STEP 3: DERIVE vs ANCESTRAL STATE CLASSIFICATION")
print("="*80)

# For each parsimony-informative site, determine state for each hominin
site_classifications = []

for pos, chars in pars_info_sites:
    ancestral = ancestral_states[pos]
    derived_states = set(chars[taxon] for taxon in ['Sima_430ka', 'Denisova3_70ka', 
                                                      'Neandertal_Vindija_44ka', 
                                                      'Feldhofer1_40ka', 
                                                      'Mezmaiskaya1_60ka'] 
                        if chars[taxon] != ancestral)
    
    if derived_states:  # Only include if there are actual derived states
        site_info = {
            'position': pos,
            'ancestral': ancestral,
            'Pan_outgroup': chars['Pan_outgroup'],
            'Sima_430ka': chars['Sima_430ka'],
            'Denisova3_70ka': chars['Denisova3_70ka'],
            'Neandertal_Vindija_44ka': chars['Neandertal_Vindija_44ka'],
            'Feldhofer1_40ka': chars['Feldhofer1_40ka'],
            'Mezmaiskaya1_60ka': chars['Mezmaiskaya1_60ka'],
            'derived_states': derived_states,
        }
        site_classifications.append(site_info)

print(f"Sites with derived states: {len(site_classifications)}")

# ============================================================================
# STEP 4: Count state changes under H1 vs H2
# ============================================================================
print("\n" + "="*80)
print("STEP 4: STATE CHANGE COUNTING (Character Parsimony)")
print("="*80)

def count_changes_h1(site_info):
    """
    H1 topology: ((Sima, (Denisova, (Vindija, Feldhofer, Mezmaiskaya))), Pan)
    Count minimum state changes needed
    """
    ancestral = site_info['ancestral']
    sima = site_info['Sima_430ka']
    denisova = site_info['Denisova3_70ka']
    vindija = site_info['Neandertal_Vindija_44ka']
    feldhofer = site_info['Feldhofer1_40ka']
    mezmaiskaya = site_info['Mezmaiskaya1_60ka']
    
    # Start with Pan (ancestral)
    # Count changes in H1 topology
    
    changes = 0
    
    # Change from ancestral to Sima lineage?
    if sima != ancestral:
        changes += 1
    
    # Denisova lineage
    if denisova != ancestral:
        # Did change occur before Denisova split?
        if sima == denisova:
            # No change needed if Sima already derived
            pass
        else:
            changes += 1
    
    # Neandertal clade
    neandertal_states = [vindija, feldhofer, mezmaiskaya]
    neandertal_ancestral = ancestral
    
    # If any neandertals are derived from ancestral
    if any(n != ancestral for n in neandertal_states):
        if denisova == ancestral:
            # Need change in Neandertal clade
            changes += 1
        else:
            # Check if change already accounted for
            if sima == denisova and denisova != ancestral:
                # No additional change needed if already derived
                pass
            else:
                changes += 1
    
    # Within Neandertals (if they differ from each other)
    if len(set(neandertal_states)) > 1:
        # Additional changes within neandertal clade
        changes += len(set(neandertal_states)) - 1
    
    return changes

def count_changes_h2(site_info):
    """
    H2 topology: (Sima, (Denisova, (Vindija, Feldhofer, Mezmaiskaya)))
    Sima in ancestral position
    Count minimum state changes needed
    """
    ancestral = site_info['ancestral']
    sima = site_info['Sima_430ka']
    denisova = site_info['Denisova3_70ka']
    vindija = site_info['Neandertal_Vindija_44ka']
    feldhofer = site_info['Feldhofer1_40ka']
    mezmaiskaya = site_info['Mezmaiskaya1_60ka']
    
    changes = 0
    
    # If Sima differs from ancestral, it's a derived state
    if sima != ancestral:
        changes += 1
    
    # Denisova branch
    if denisova != ancestral:
        if denisova != sima:
            changes += 1
    
    # Neandertal clade
    neandertal_states = [vindija, feldhofer, mezmaiskaya]
    
    if any(n != ancestral for n in neandertal_states):
        if all(n == denisova for n in neandertal_states):
            # No change if same as Denisova
            pass
        else:
            changes += 1
    
    # Within Neandertals
    if len(set(neandertal_states)) > 1:
        changes += len(set(neandertal_states)) - 1
    
    return changes

# Count changes for all sites
h1_total = 0
h2_total = 0
h1_better = 0
h2_better = 0
h1_equal = 0

for site in site_classifications:
    h1 = count_changes_h1(site)
    h2 = count_changes_h2(site)
    
    h1_total += h1
    h2_total += h2
    
    if h1 < h2:
        h1_better += 1
    elif h2 < h1:
        h2_better += 1
    else:
        h1_equal += 1

print(f"\nCharacter-based parsimony counts:")
print(f"  H1 (Sima sister-group):  {h1_total} total state changes")
print(f"  H2 (Sima ancestral):     {h2_total} total state changes")
print(f"  Difference:              {abs(h1_total - h2_total)} steps")

print(f"\nSite-by-site comparison:")
print(f"  H1 better:  {h1_better} sites")
print(f"  H2 better:  {h2_better} sites")
print(f"  Tied:       {h1_equal} sites")

print("\n" + "="*80)
if h1_total < h2_total:
    improvement = 100 * (h2_total - h1_total) / h2_total
    print(f"✓ H1 (sister-group) is MORE PARSIMONIOUS")
    print(f"  Reduction: {improvement:.1f}%")
elif h2_total < h1_total:
    improvement = 100 * (h1_total - h2_total) / h1_total
    print(f"✓ H2 (ancestral) is MORE PARSIMONIOUS")
    print(f"  Reduction: {improvement:.1f}%")
else:
    print(f"= Both topologies require equal state changes (tied)")
print("="*80)

# ============================================================================
# STEP 5: Validation with sample sites
# ============================================================================
print("\n" + "="*80)
print("STEP 5: EXAMPLE SITES (first 10 parsimony-informative)")
print("="*80)

for i, site in enumerate(site_classifications[:10]):
    print(f"\nSite {i+1} (position {site['position']}):")
    print(f"  Ancestral state (from Pan): {site['ancestral']}")
    print(f"  Sima:            {site['Sima_430ka']}")
    print(f"  Denisova3:       {site['Denisova3_70ka']}")
    print(f"  Vindija:         {site['Neandertal_Vindija_44ka']}")
    print(f"  Feldhofer:       {site['Feldhofer1_40ka']}")
    print(f"  Mezmaiskaya:     {site['Mezmaiskaya1_60ka']}")
    
    h1 = count_changes_h1(site)
    h2 = count_changes_h2(site)
    winner = "H1" if h1 < h2 else ("H2" if h2 < h1 else "TIE")
    print(f"  H1 steps: {h1}  |  H2 steps: {h2}  →  {winner}")

print("\n" + "="*80)
print("ANALYSIS COMPLETE")
print("="*80)
