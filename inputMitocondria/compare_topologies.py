#!/usr/bin/env python3
################################################################################
# COMPARE_TOPOLOGIES.py
# 
# Compares tree topologies from IQ-TREE, RAxML, and MrBayes
# Extracts key clades and calculates topology concordance
#
# Usage: python3 compare_topologies.py
# Run from: ~/sima/inputMitocondria/phase3a_results/
################################################################################

import re
from pathlib import Path
from collections import defaultdict

def extract_tree_string(filepath):
    """Extract Newick tree string from output file"""
    try:
        with open(filepath) as f:
            content = f.read()
        
        # Find tree string (last line with Newick format)
        lines = content.strip().split('\n')
        for line in reversed(lines):
            if line.startswith('(') and line.endswith(';'):
                return line
        return None
    except:
        return None

def extract_clade(tree_str, taxon):
    """Find the clade position of a taxon in tree"""
    try:
        # Simple clade extraction - find parent positions
        idx = tree_str.find(taxon)
        if idx == -1:
            return None
        
        # Count parentheses to find clade depth
        depth = 0
        for i in range(idx):
            if tree_str[i] == '(':
                depth += 1
            elif tree_str[i] == ')':
                depth -= 1
        return depth
    except:
        return None

def extract_bipartitions(tree_str):
    """Extract bipartitions (major splits) from tree"""
    bipartitions = []
    
    # Simple regex to find major groups
    # Look for patterns like (A,B) indicating sister taxa
    taxa_patterns = [
        ('Sima', 'Neandertal'),
        ('Sima', 'Denisova'),
        ('Neandertal', 'Denisova'),
        ('Feldhofer', 'Mezmaiskaya'),
    ]
    
    result = {}
    for tax1, tax2 in taxa_patterns:
        # Check if they're in same parentheses
        pattern = f'\\(.*{tax1}.*{tax2}.*\\)|\\(.*{tax2}.*{tax1}.*\\)'
        if re.search(pattern, tree_str):
            result[f"{tax1}+{tax2}"] = True
        else:
            result[f"{tax1}+{tax2}"] = False
    
    return result

def extract_support_values(tree_str):
    """Extract support values (bootstrap/posterior) from tree"""
    # Find all values in brackets or after colons
    support_vals = {
        'bootstrap': re.findall(r'\):(\d+)', tree_str),
        'posterior': re.findall(r'\[&prob=([\d.]+)\]', tree_str),
    }
    return support_vals

def main():
    workdir = Path('.')
    
    print("=" * 120)
    print("PHASE 3a TOPOLOGY COMPARISON")
    print("=" * 120)
    print()
    
    # Collect tree files
    trees = {
        'IQ-TREE (HKY+G4)': workdir / "iqtree_hky.treefile",
        'RAxML-NG (GTR+G)': workdir / "raxml_gtr.raxml.bestTree",
        'MrBayes (HKY+G)': workdir / "mrbayes.con.tre",
    }
    
    tree_data = {}
    
    for method, filepath in trees.items():
        print(f"{method}")
        print("-" * 120)
        
        if not filepath.exists():
            print(f"  ⏳ Results not ready yet (run phase3a_master.sh first)")
            print()
            continue
        
        tree_str = extract_tree_string(str(filepath))
        if not tree_str:
            print(f"  ⚠️  Could not parse tree from {filepath.name}")
            print()
            continue
        
        tree_data[method] = tree_str
        
        print(f"  Tree: {tree_str[:100]}...")
        print()
        
        # Extract key bipartitions
        biparts = extract_bipartitions(tree_str)
        print("  Key Bipartitions (Sister Taxa):")
        for pair, grouped in biparts.items():
            status = "✓ Grouped" if grouped else "✗ Separate"
            print(f"    {pair:.<30} {status}")
        
        print()
        
        # Extract support values
        supports = extract_support_values(tree_str)
        if supports['bootstrap']:
            bs_vals = [int(x) for x in supports['bootstrap']]
            print(f"  Bootstrap Support:")
            print(f"    Min:   {min(bs_vals)}")
            print(f"    Max:   {max(bs_vals)}")
            print(f"    Mean:  {sum(bs_vals)/len(bs_vals):.1f}")
            print(f"    ≥95%:  {sum(1 for x in bs_vals if x >= 95)}/{len(bs_vals)}")
        
        if supports['posterior']:
            pp_vals = [float(x) for x in supports['posterior']]
            print(f"  Posterior Probability:")
            print(f"    Min:   {min(pp_vals):.3f}")
            print(f"    Max:   {max(pp_vals):.3f}")
            print(f"    Mean:  {sum(pp_vals)/len(pp_vals):.3f}")
            print(f"    ≥0.95: {sum(1 for x in pp_vals if x >= 0.95)}/{len(pp_vals)}")
        
        print()
    
    # ========== TOPOLOGY CONCORDANCE ==========
    print("=" * 120)
    print("TOPOLOGY CONCORDANCE ANALYSIS")
    print("=" * 120)
    print()
    
    if len(tree_data) >= 2:
        print("Comparing topologies across methods:")
        print()
        
        # Extract key clade positions
        key_splits = {
            'Sima basal position': None,
            'Neandertal+Denisova clade': None,
            'Feldhofer+Mezmaiskaya clade': None,
        }
        
        for method, tree_str in tree_data.items():
            print(f"  {method}:")
            
            # Check major splits
            if re.search(r'\(Sima.*,\(.*Neandertal.*Denisova', tree_str):
                print("    → Sima BASAL to (Neandertal + Denisova) [mtDNA H0 topology]")
            elif re.search(r'\(\(Sima.*Neandertal.*\).*Denisova', tree_str):
                print("    → (Sima + Neandertal) sister to Denisova [nDNA H1 topology]")
            else:
                print("    → Other topology")
            
            print()
    else:
        print("  ⏳ Need at least 2 methods completed for comparison")
        print()
    
    # ========== KH TEST SUMMARY ==========
    print("=" * 120)
    print("TOPOLOGY HYPOTHESIS TESTS (KH/AU)")
    print("=" * 120)
    print()
    
    h0_file = workdir / "khtest_h0_mtdna.iqtree"
    h1_file = workdir / "khtest_h1_ndna.iqtree"
    
    hypothesis_results = {}
    
    for test_name, filepath in [("H0 (mtDNA)", h0_file), ("H1 (nDNA)", h1_file)]:
        if filepath.exists():
            try:
                with open(filepath) as f:
                    content = f.read()
                
                print(f"{test_name} - {filepath.name}")
                print("-" * 120)
                
                # Extract AU test p-value
                au_match = re.search(r'AU test\s+p-value\s*=\s+([\d.]+)', content)
                if au_match:
                    au_pval = float(au_match.group(1))
                    hypothesis_results[test_name] = au_pval
                    status = "NOT REJECTED" if au_pval > 0.05 else "REJECTED"
                    print(f"  AU test p-value: {au_pval:.4f} → {status}")
                
                # Extract KH test p-value
                kh_match = re.search(r'KH test\s+p-value\s*=\s+([\d.]+)', content)
                if kh_match:
                    kh_pval = float(kh_match.group(1))
                    status = "NOT REJECTED" if kh_pval > 0.05 else "REJECTED"
                    print(f"  KH test p-value: {kh_pval:.4f} → {status}")
                
                # Extract tree likelihood
                ll_match = re.search(r'Log-likelihood.*?:\s+([-\d.]+)', content)
                if ll_match:
                    print(f"  Log-likelihood: {ll_match.group(1)}")
                
                print()
            except Exception as e:
                print(f"  Error reading {filepath.name}: {e}")
                print()
        else:
            print(f"{test_name}: ⏳ Results not ready yet")
            print()
    
    if hypothesis_results:
        print()
        print("HYPOTHESIS COMPARISON:")
        print("-" * 120)
        
        h0_pval = hypothesis_results.get("H0 (mtDNA)")
        h1_pval = hypothesis_results.get("H1 (nDNA)")
        
        if h0_pval and h1_pval:
            if h0_pval > h1_pval:
                print(f"  H0 (mtDNA) is BETTER SUPPORTED (p={h0_pval:.4f} > p={h1_pval:.4f})")
            else:
                print(f"  H1 (nDNA) is BETTER SUPPORTED (p={h1_pval:.4f} > p={h0_pval:.4f})")
            
            if abs(h0_pval - h1_pval) < 0.01:
                print(f"  ⚠️  Topologies are statistically similar (Δp < 0.01)")
            else:
                print(f"  ✓ Clear preference for one topology (Δp = {abs(h0_pval - h1_pval):.4f})")
    
    print()
    print("=" * 120)
    print()

if __name__ == "__main__":
    main()
