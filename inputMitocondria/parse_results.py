#!/usr/bin/env python3
################################################################################
# PARSE_RESULTS.py - Extract bootstrap/posterior values from Phase 3a outputs
# 
# Usage: python3 parse_results.py
# Run from: ~/sima/inputMitocondria/phase3a_results/
#
# Parses:
#   - IQ-TREE .iqtree file (bootstrap support)
#   - RAxML-NG .raxml.support file (bootstrap support)
#   - MrBayes .con.tre file (posterior probabilities)
#   - KH test .iqtree files (AU p-values)
################################################################################

import re
import os
from pathlib import Path

def parse_newick_support(treefile):
    """Extract bootstrap/posterior values from Newick tree string"""
    try:
        with open(treefile) as f:
            tree = f.read().strip()
        
        # Find all support values in brackets
        supports = re.findall(r'\):([0-9.]+)', tree)
        
        if supports:
            # Convert to percentages if needed
            nums = [float(x) for x in supports]
            if max(nums) <= 1:
                nums = [x * 100 for x in nums]
            return nums
        return None
    except:
        return None

def parse_iqtree_file(filepath):
    """Extract key statistics from IQ-TREE output"""
    try:
        with open(filepath) as f:
            content = f.read()
        
        result = {}
        
        # Log-likelihood
        ll_match = re.search(r'Log-likelihood of the tree:\s+([-\d.]+)', content)
        if ll_match:
            result['log_likelihood'] = ll_match.group(1)
        
        # Model
        model_match = re.search(r'Model of substitution:\s+(\S+)', content)
        if model_match:
            result['model'] = model_match.group(1)
        
        # Bootstrap count
        if '1000 bootstrap trees' in content:
            result['bootstrap'] = 1000
        
        # Tree topology (extract key clade info)
        tree_match = re.search(r'(.*);', content)
        if tree_match:
            tree = tree_match.group(1)
            # Check for Sima position
            if 'Sima' in tree:
                if tree.find('Sima') < tree.find('Neandertal'):
                    result['sima_position'] = 'basal to Neandertal'
                else:
                    result['sima_position'] = 'grouped with Neandertal'
        
        return result
    except Exception as e:
        print(f"Error parsing {filepath}: {e}")
        return {}

def parse_raxml_file(filepath):
    """Extract key statistics from RAxML output"""
    try:
        with open(filepath) as f:
            content = f.read()
        
        result = {}
        
        # Final log-likelihood
        ll_match = re.search(r'Final GAMMA-based Score of best tree\s*:\s+([-\d.]+)', content)
        if ll_match:
            result['log_likelihood'] = ll_match.group(1)
        
        # Bootstrap info
        if 'bootstrap' in content.lower():
            result['bootstrap'] = 1000
        
        return result
    except:
        return {}

def parse_mrbayes_tree(treefile):
    """Extract posterior probabilities from MrBayes consensus tree"""
    try:
        with open(treefile) as f:
            tree = f.read().strip()
        
        # Extract posterior values
        posteriors = re.findall(r'\[&prob=([\d.]+)', tree)
        if posteriors:
            posteriors = [float(x) for x in posteriors]
            return posteriors
        return None
    except:
        return None

def parse_khtest_file(filepath):
    """Extract AU test p-values from KH test output"""
    try:
        with open(filepath) as f:
            content = f.read()
        
        result = {}
        
        # AU test p-values
        au_match = re.search(r'AU test\s+p-value\s*=\s+([\d.]+)', content)
        if au_match:
            result['au_pvalue'] = au_match.group(1)
        
        # KH test p-values
        kh_match = re.search(r'KH test\s+p-value\s*=\s+([\d.]+)', content)
        if kh_match:
            result['kh_pvalue'] = kh_match.group(1)
        
        # Tree likelihood
        ll_match = re.search(r'Log-likelihood.*?:\s+([-\d.]+)', content)
        if ll_match:
            result['log_likelihood'] = ll_match.group(1)
        
        return result
    except Exception as e:
        print(f"Error parsing {filepath}: {e}")
        return {}

def main():
    workdir = Path('.')
    
    print("=" * 100)
    print("PHASE 3a RESULTS PARSING & SUMMARY")
    print("=" * 100)
    print()
    
    # ========== IQ-TREE RESULTS ==========
    print("IQ-TREE (HKY+G4) - Maximum Likelihood + Bootstrap")
    print("-" * 100)
    
    iqtree_file = workdir / "iqtree_hky.iqtree"
    if iqtree_file.exists():
        iqtree_results = parse_iqtree_file(str(iqtree_file))
        if iqtree_results:
            for key, val in iqtree_results.items():
                print(f"  {key:.<30} {val}")
            
            # Extract bootstrap support values
            tree_file = workdir / "iqtree_hky.treefile"
            if tree_file.exists():
                supports = parse_newick_support(str(tree_file))
                if supports:
                    print(f"  {'bootstrap_values':.<30} min={min(supports):.1f}, max={max(supports):.1f}, mean={sum(supports)/len(supports):.1f}")
                    print(f"  {'nodes_with_100_support':.<30} {sum(1 for s in supports if s >= 100)}/{len(supports)}")
    else:
        print("  ⏳ Waiting for results (run phase3a_master.sh first)")
    
    print()
    
    # ========== RAxML-NG RESULTS ==========
    print("RAxML-NG (GTR+G) - Maximum Likelihood + Bootstrap")
    print("-" * 100)
    
    raxml_file = workdir / "raxml_gtr.raxml.log"
    if raxml_file.exists():
        raxml_results = parse_raxml_file(str(raxml_file))
        if raxml_results:
            for key, val in raxml_results.items():
                print(f"  {key:.<30} {val}")
            
            # Extract bootstrap support
            support_file = workdir / "raxml_gtr.raxml.support"
            if support_file.exists():
                supports = parse_newick_support(str(support_file))
                if supports:
                    print(f"  {'bootstrap_values':.<30} min={min(supports):.1f}, max={max(supports):.1f}, mean={sum(supports)/len(supports):.1f}")
                    print(f"  {'nodes_with_100_support':.<30} {sum(1 for s in supports if s >= 100)}/{len(supports)}")
    else:
        print("  ⏳ Waiting for results (run phase3a_master.sh first)")
    
    print()
    
    # ========== MrBayes RESULTS ==========
    print("MrBayes (HKY+G) - Bayesian MCMC Inference")
    print("-" * 100)
    
    mrbayes_tree = workdir / "mrbayes.con.tre"
    if mrbayes_tree.exists():
        posteriors = parse_mrbayes_tree(str(mrbayes_tree))
        if posteriors:
            print(f"  {'posterior_probabilities':.<30} min={min(posteriors):.3f}, max={max(posteriors):.3f}, mean={sum(posteriors)/len(posteriors):.3f}")
            print(f"  {'nodes_with_pp_>0.95':.<30} {sum(1 for p in posteriors if p >= 0.95)}/{len(posteriors)}")
        
        # Extract summary statistics
        mrbayes_sum = workdir / "mrbayes.sum"
        if mrbayes_sum.exists():
            try:
                with open(mrbayes_sum) as f:
                    content = f.read()
                    if "MCMC chain" in content:
                        print(f"  {'MCMC_status':.<30} Completed")
            except:
                pass
    else:
        print("  ⏳ Waiting for results (run phase3a_master.sh first)")
    
    print()
    
    # ========== KH TOPOLOGY TESTS ==========
    print("KISHINO-HASEGAWA (KH) TOPOLOGY TESTS")
    print("-" * 100)
    
    # H0 test (mtDNA topology)
    h0_file = workdir / "khtest_h0_mtdna.iqtree"
    if h0_file.exists():
        h0_results = parse_khtest_file(str(h0_file))
        print(f"H0 (mtDNA topology): Sima basal to (Denisova + Neandertal)")
        for key, val in h0_results.items():
            print(f"  {key:.<30} {val}")
    else:
        print(f"H0 results: ⏳ Waiting for results")
    
    print()
    
    # H1 test (nDNA topology)
    h1_file = workdir / "khtest_h1_ndna.iqtree"
    if h1_file.exists():
        h1_results = parse_khtest_file(str(h1_file))
        print(f"H1 (nDNA topology): (Sima + Neandertal) sister to Denisova")
        for key, val in h1_results.items():
            print(f"  {key:.<30} {val}")
    else:
        print(f"H1 results: ⏳ Waiting for results")
    
    print()
    print("=" * 100)
    print()
    
    # ========== INTERPRETATION ==========
    print("INTERPRETATION:")
    print("-" * 100)
    print()
    print("Bootstrap/Posterior Support Thresholds:")
    print("  ≥95%  → Strong support (confident clade)")
    print("  80-94% → Moderate support (well-supported)")
    print("  <80%  → Weak support (uncertain clade)")
    print()
    print("Topology Hypothesis Tests:")
    print("  AU p-value > 0.05 → Topology not significantly rejected")
    print("  AU p-value ≤ 0.05 → Topology significantly rejected")
    print()
    print("Expected Results:")
    print("  - Both IQ-TREE and RAxML should show ≥95% bootstrap for main clades")
    print("  - MrBayes posterior probability ≥0.95 for well-supported clades")
    print("  - KH tests will show which topology (H0 vs H1) is better supported")
    print()
    
    print("=" * 100)
    print("Next steps: Review trees, extract specific support values for manuscript")
    print("=" * 100)

if __name__ == "__main__":
    main()
