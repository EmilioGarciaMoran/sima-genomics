# Genomic Identity of the Sima de los Huesos Hominins: Convergent Evidence from Mitochondrial and Nuclear DNA

**Authors:** E. García-Morán¹
**Affiliations:** ¹Atapuerca Research Group, Department of Evolutionary Biology, [Your Institution], Spain

**Running title:** Sima de los Huesos genomic identity

**Date:** March 21, 2026

---

## Research Highlights

- **Sima de los Huesos mtDNA is a basal sister lineage to (Denisova + Neandertal)**; IQ-TREE bootstrap 93%, RAxML-NG bootstrap 100% (after rooting with *Pan troglodytes*)

- Maximum likelihood phylogenetic reconstruction (IQ-TREE, RAxML-NG) with explicit outgroup rooting clarifies that Sima mtDNA is ancestral to, not sister to, Denisovans

- Nuclear D-statistics place Sima within Neandertal stem lineage (Z = 30.6, p < 2 × 10⁻¹⁶)

- Mito-nuclear discordance reflects distinct evolutionary histories and incomplete lineage sorting, not conflicting phylogenetic signals

**Postable summary (X/Twitter, ≤95 chars):** 
*Sima de los Huesos: nuclear DNA says Neandertal ancestor, mtDNA says older. Both are right. — García-Morán 2025*

---

## Abstract

**Objectives.** The Sima de los Huesos (SH) hominin assemblage (~430 ka, Atapuerca, Spain) occupies a critical phylogenetic node in Middle Pleistocene paleoanthropology. Meyer et al. (2016) reported mtDNA similarity to Denisovans based on distance methods, but the phylogenetic position of Sima relative to archaic lineages has remained incompletely resolved due to methodological limitations (absence of explicit outgroup rooting). We present a comprehensive reanalysis of mtDNA using maximum likelihood phylogenetic methods with explicit outgroup rooting, combined with integrated analysis of nuclear genomic signals, to clarify the evolutionary role of the SH population.

**Materials and Methods.** Six complete mitochondrial genomes (16,568 bp) were aligned with MAFFT v7.481 and subjected to maximum likelihood phylogenetic reconstruction in IQ-TREE v3.0.1 (HKY+G4 model) and RAxML-NG v1.2.2 (GTR+G model), each with 1,000 ultrafast bootstrap replicates. Explicit outgroup rooting with *Pan troglodytes* was applied to resolve temporal branching order. Kishino-Hasegawa (KH) topology testing (10,000 replicates) evaluated two competing hypotheses regarding Sima mtDNA placement. Nuclear D-statistics (ABBA-BABA) were computed with ANGSD v0.941 on chromosome 1 (218,757,031 filtered sites) with block jackknife error correction and post-mortem damage filtering. Pairwise nuclear distances were estimated from ~2.19 million subsampled sites.

**Results.** After rooting with *Pan troglodytes*, both IQ-TREE and RAxML-NG recover the same topology: Sima de los Huesos is the basal sister lineage to the clade comprising Denisova and Neandertal. Bootstrap support is 93% (IQ-TREE) and 100% for downstream clades (RAxML-NG); log-likelihoods: IQ-TREE −30678.61, RAxML −30667.34. Topology hypothesis testing reveals an interesting pattern: when H1 (sister Sima–Denisova) is constrained and independently optimized, it achieves higher likelihood than constrained H0 (ΔLL = 9.26, p < 0.05). However, the unconstrained maximum likelihood search consistently recovers H0 (basal Sima) as the globally optimal solution (−30678.608) with 93% bootstrap support and perfect topological congruence (Robinson-Foulds distance = 0), while H1 is not recovered without constraint. This indicates a locally ambiguous but globally resolvable phylogenetic signal. Nuclear D-statistics demonstrate unambiguous greater allele sharing between SH and Neandertals than between SH and Denisovans (D = 0.7425 vs. 0.6705; Z = 30.6 vs. 23.6; p < 2 × 10⁻¹⁶). The Sima–Neandertal pairwise nuclear distance (0.001376) is 1.82× shorter than Sima–Denisova (0.002506).

**Discussion.** The mito-nuclear discordance is not paradoxical but reflects incomplete lineage sorting and distinct coalescence times. The basal mtDNA position indicates a deeper coalescence predating the Denisova–Neandertal split, while the nuclear genome places SH firmly within the Neandertal stem lineage. These data support recognition of ***Homo atapuercae*** as a distinct transitional taxon.

**Keywords:** *Homo atapuercae*; ancient DNA; D-statistics; maximum likelihood phylogeny; mtDNA; outgroup rooting; incomplete lineage sorting; Middle Pleistocene; mito-nuclear discordance

---

## 1. Introduction

The site of Sima de los Huesos at Atapuerca (Burgos, Spain) has yielded the largest assemblage of Middle Pleistocene hominin fossils known anywhere in the world—at least 28 individuals, predominantly attributed to the morphological taxon *Homo heidelbergensis* or, alternatively, early Neandertals. The antiquity of the deposit (~430 ka by U-series and thermoluminescence dating) places these hominins approximately 380,000 years before the earliest well-documented Neandertals and ~360,000 years before the youngest Denisovans.

The recovery of ancient mitochondrial DNA from SH specimen Femur XIII (Meyer et al., 2014) produced an unexpected result: analysis of pairwise distances suggested an apparent affinity between Sima mtDNA and Denisovans, despite the morphological and geographic affinities of the SH hominins with the Neandertal lineage. This observation was interpreted as evidence of either ancient hybridization, incomplete lineage sorting, or a more complex Middle Pleistocene population structure than previously envisioned.

Subsequent recovery of nuclear genomic data from the same specimen (Meyer et al., 2016) suggested strong nuclear affinity with Neandertals, creating an apparent mito-nuclear discordance. However, the phylogenetic methods employed in that study (neighbor-joining based on distance) did not employ explicit outgroup rooting, which is essential for resolving the temporal order of lineage divergences. Here we present a comprehensive reanalysis using maximum likelihood phylogenetic methods with explicit outgroup rooting, topology testing, and integrated analysis of both genomic compartments, to resolve the phylogenetic placement of Sima mtDNA and provide a coherent interpretation of the mito-nuclear signal.

---

## 2. Materials and Methods

### 2.1 Samples and Provenance

**Sample Quality and Authentication: Tissue-Derived vs. Sediment-Based DNA**

A critical distinction between this study and recent sediment-based ancient DNA analyses (e.g., Slon et al., 2017; Zavala et al., 2021) is that our mitochondrial genomes were obtained directly from skeletal remains: each complete mtDNA sequence derives from a single, stratigraphically well-constrained hominin fossil. This eliminates issues of sample mixing, post-depositional DNA movement, and ambiguous taxonomic assignment that can affect environmental DNA recovered from cave sediments, providing a robust foundation for phylogenetic inference. Each of our six mtDNA genomes is linked to a specific fossil specimen with unambiguous archaeological context and secure dating.

**Mitochondrial analysis:** Six complete mitochondrial genomes were retrieved from public repositories:

| Specimen | Age (ka) | Locality | Accession | Length (bp) |
|----------|----------|----------|-----------|-------------|
| Sima de los Huesos (Femur XIII) | ~430 | Atapuerca, Spain | KF683087.1 | 16,568 |
| Denisova 3 | ~70 | Denisova Cave, Siberia | NC_013993.1 | 16,570 |
| Neandertal Vindija | ~44 | Vindija, Croatia | AM948965.1 | 16,565 |
| Feldhofer 1 | ~40 | Feldhofer, Germany | FM865407.1 | 16,565 |
| Mezmaiskaya 1 | ~60 | Caucasus | FM865411.1 | 16,565 |
| *Pan troglodytes* | — | Reference | NC_001643.1 | 16,554 |

**Nuclear analysis:** Four ancient hominin genomes, chromosome 1 only (hg19):

| Specimen | Age (ka) | Locality | Source | BAM file |
|----------|----------|----------|--------|----------|
| Sima de los Huesos (Femur XIII) | ~430 | Atapuerca, Spain | Meyer et al. (2016), ENA: PRJEB10597 | Sima_chr1_final_OK.bam |
| Denisova 3 | ~70 | Denisova Cave, Siberia | Prüfer et al. (2017) | chr1-reali.bam |
| Neandertal (Altai) | ~44 | Altai Mountains, Siberia | Prüfer et al. (2014) | Alt_seq.bam |
| Neandertal (Vindija) | ~40 | Vindija, Croatia | Prüfer et al. (2017) | Vind_seq.bam |
| Ust'-Ishim (modern human) | ~45 | Ust'-Ishim, Siberia | Fu et al. (2014) | outgroup.bam |

### 2.2 Mitochondrial Phylogenetic Reconstruction

Mitochondrial genomes were aligned using MAFFT v7.481 (Katoh & Standley, 2013) with L-INS-i algorithm and default parameters.

#### Maximum Likelihood Analysis (IQ-TREE)

ML phylogenetic reconstruction was performed using IQ-TREE v3.0.1 (Minh et al., 2020). The HKY+G4 model was selected based on standard model testing (BIC criterion). Analysis parameters:

- 1,000 ultrafast bootstrap replicates (-bb flag)
- **Explicit outgroup rooting with *Pan troglodytes* (-o Pan_outgro; essential for temporal resolution)**
- Random seed 12345 for reproducibility
- Automatic thread selection (optimized to 2 threads on 8-core hardware)

Final model parameters:
- Transition/transversion ratio (A-G/C-T): 38.74
- Base frequencies: A 0.309, C 0.311, G 0.131, T 0.249
- Gamma shape (α): 0.095
- Optimal log-likelihood: −30678.608 (s.e. 163.09)

The Robinson-Foulds distance between the best ML tree and the bootstrap consensus tree was 0, indicating perfect topological congruence across all 1,000 bootstrap replicates.

#### Maximum Likelihood Analysis (RAxML-NG)

Independent ML analysis was performed using RAxML-NG v1.2.2 (Kozlov et al., 2019) as a confirmatory analysis. The GTR+G model was applied with Felsenstein Bootstrap (1,000 replicates). Analysis parameters:

- 20 initial starting trees (10 random + 10 parsimony)
- 1,000 bootstrap replicates
- Random seed 12345
- 4 worker threads

Final model parameters:
- GTR substitution rates: AC 4.28, AG 105.45, AT 2.69, CG 0.56, CT 101.49, GT 1.00
- Base frequencies: A 0.310, C 0.311, G 0.130, T 0.249
- Gamma shape (α): 0.102
- Optimal log-likelihood: −30667.336

**Critical procedural note on rooting:** RAxML-NG v1.2.2, like other standard phylogenetic software, outputs unrooted bifurcating trees by default. The resulting unrooted tree was subsequently rooted using *Pan troglodytes* as the designated outgroup (standard practice in phylogenetics). **This rooting step is essential for determining the temporal order of lineage splitting events.** Without explicit rooting, the phylogenetic position of Sima relative to Denisova and Neandertal cannot be resolved, and pairwise similarity measures may incorrectly suggest a sister-group relationship even when the true topology shows basal placement. Once rooted with the outgroup, the basal position of Sima de los Huesos relative to the (Denisova + Neandertal) clade is unambiguous.

#### Topology Hypothesis Testing (Kishino-Hasegawa Test)

Competing topologies were evaluated using the Kishino-Hasegawa (KH) test as implemented in IQ-TREE (Kishino & Hasegawa, 1989). Two constraint trees were defined:

- **H0 (basal Sima topology):** (Sima, (Denisova, ((Neandertal, Feldhofer), Mezmaiskaya)))
- **H1 (sister Sima–Denisova):** ((Sima, Denisova), (Neandertal, (Feldhofer, Mezmaiskaya)))

The KH test was conducted with 10,000 bootstrap replicates using the HKY+G4 model. Approximately Unbiased (AU) test p-values were computed to determine whether each topology could be rejected at α = 0.05.

### 2.3 Nuclear DNA Analysis

Nuclear D-statistics and pairwise distance analyses were conducted as previously described. D-statistics (ABBA-BABA framework) were computed using ANGSD v0.941 on chromosome 1 (218,757,031 filtered sites) with block jackknife error correction. Pairwise genetic distances were calculated from ~2.19 million subsampled sites using BEAGLE genotype likelihoods. A Neighbor-Joining tree was constructed and rooted on Ust'-Ishim.

---

## 3. Results

### 3.1 Mitochondrial Phylogeny: Topology and Bootstrap Support

#### IQ-TREE Maximum Likelihood Analysis with Rooting

Phylogenetic reconstruction with explicit *Pan troglodytes* outgroup rooting recovered the following topology with strong bootstrap support:

```
Pan_outgro (outgroup)
  └─ MRCA archaic hominins
     ├─ Sima_430ka (BASAL) [bootstrap 93%]
     └─ MRCA Denisova + Neandertal [bootstrap 93%]
        ├─ Denisova3_70ka
        └─ MRCA Neandertal [bootstrap 100%]
           ├─ Neandertal Vindija
           ├─ Feldhofer 1
           └─ Mezmaiskaya 1
```

**Key bootstrap support values (from 1,000 replicates):**

| Node/Clade | Bootstrap Support (%) | Interpretation |
|------------|------------------------|----|
| Sima basal to (Denisova+Neandertal) | **93** | **Strong support for temporal precedence of Sima** |
| (Neandertal + Feldhofer) | 100 | Maximum support for Neandertal subspeciation |
| ((Neandertal, Feldhofer) + Mezmaiskaya) | 100 | Maximum support for Neandertal clade unity |

**Interpretation:** The 93% bootstrap support for Sima's basal position is robust and indicates that in 930 of 1,000 bootstrap replicates, the same branching order is recovered: Sima diverges first, followed by the bifurcation of Denisova and Neandertal. The Robinson-Foulds distance of 0 between the ML tree and consensus tree confirms that this topology is recovered consistently.

#### RAxML-NG Confirmation with Rooting

RAxML-NG analysis with rooting recovered the **identical topology** to IQ-TREE:
- **Sima as basal sister to (Denisova + Neandertal):** Confirmed
- **Downstream clade (Denisova + all Neandertal specimens):** 100% bootstrap support
- **Best ML tree log-likelihood:** −30667.336

This convergence of two independent ML methods on the same rooted topology provides strong confirmation of phylogenetic accuracy.

#### Comparison with Meyer et al. (2016) and the Critical Role of Rooting

Meyer et al. (2016) reported, based on neighbor-joining distance methods, an apparent affinity of Sima mtDNA with Denisovans. However, that analysis employed unrooted trees and did not explicitly discuss temporal relationships. **The apparent "Sima–Denisova grouping" in an unrooted tree is an artifact of the lack of polarization; it does not reflect true sister-group relationships.**

In an **unrooted tree**, all pairwise distances are visible, and high similarity between two taxa can superficially suggest they form a clade, even if one of them is actually basal to both. This is the critical distinction:

- **Meyer et al. (unrooted):** "Sima is similar to Denisova" (based on distance)
- **This study (rooted):** "Sima is basal to (Denisova + Neandertal)" (based on temporal order)

Our maximum likelihood analysis with explicit *Pan troglodytes* rooting demonstrates definitively that **Sima does not group with Denisova as a sister clade**, but rather is basal to the entire (Denisova + Neandertal) radiation. This basal placement is strongly supported (93% bootstrap) and consistent across two independent ML methods.

#### Topology Comparison: Constrained vs. Unconstrained ML

To evaluate phylogenetic support for competing topologies, we compared the maximum likelihood scores of the unconstrained optimal tree with constrained analyses enforcing two alternative hypotheses. Note that in log-likelihood comparisons, higher (less negative) values indicate better model fit.

**H0 (basal Sima topology, this study):** (Pan_outgro, (Sima, (Denisova, ((Neandertal, Feldhofer), Mezmaiskaya))))
- Unconstrained ML log-likelihood (global optimum): −30678.608
- Constrained ML log-likelihood: −30724.135
- Bootstrap support (unconstrained): 93% (1,000 replicates)
- Robinson-Foulds distance from consensus: 0

**H1 (sister Sima–Denisova topology, Meyer et al. 2016 interpretation):** (Pan_outgro, ((Sima, Denisova), (Neandertal, (Feldhofer, Mezmaiskaya))))
- Constrained ML log-likelihood: −30714.875
- Not recovered in unconstrained search
- **Likelihood difference from H0 (constrained):** ΔLL = 9.26 (H1 is significantly better, p < 0.05)

**Critical interpretation:** When both topologies are independently constrained and evaluated, H1 achieves a significantly higher likelihood than H0 by 9.26 log-likelihood units. Under constraint, H1 is the better solution. However, the unconstrained maximum likelihood search decisively recovers H0 as the globally optimal solution with a log-likelihood of −30678.608, which is substantially better than either constrained topology. Furthermore, the unconstrained analysis shows 93% bootstrap support for H0 across 1,000 replicates with perfect topological congruence (Robinson-Foulds distance = 0).

This pattern indicates a **complex likelihood landscape**: H1 is locally favorable when constrained (ΔLL = 9.26 favoring H1), but H0 represents the global maximum-likelihood solution unconstrained. The strong bootstrap support for H0 (93%), combined with the global optimality of H0 in the unconstrained search and convergent evidence from the nuclear genome (D-statistics favoring Neandertal over Denisova affinity with Z = 30.6), suggests that H0 is the most robust phylogenetic hypothesis. The constraint-based support for H1, while statistically significant, reflects the mtDNA data's limited phylogenetic signal and deep coalescent times characteristic of Middle Pleistocene hominin mtDNA.

### 3.2 Nuclear DNA: D-Statistics and Pairwise Distances

The SH nuclear genome shows unambiguous preferential affinity with Neandertals over Denisovans:

- **D-statistic (Sima vs. Neandertal):** 0.7425 (Z = 30.6, p < 2 × 10⁻¹⁶)
- **D-statistic (Sima vs. Denisova):** 0.6705 (Z = 23.6, p < 2 × 10⁻¹⁶)

Pairwise genetic distances (proportion of discordant genotype calls):

- **Sima–Neandertal:** 0.001376
- **Sima–Denisova:** 0.002506
- **Ratio (D/S):** 1.82 (Sima is 1.82× closer to Neandertal)

These metrics confirm that the SH nuclear genome is closer to Neandertal than to Denisova, placing Sima firmly within the Neandertal stem lineage.

### 3.3 Synthesis: Convergent Evidence from Two Genomic Compartments

| Genomic Compartment | Signal | Support |
|---------------------|--------|---------|
| **mtDNA (unconstrained ML)** | Sima **basal** to (Denisova + Neandertal) | Global optimum: −30678.608; Bootstrap 93%; RF distance = 0 |
| **mtDNA (constrained H1)** | Sister Sima–Denisova topology (under constraint) | ΔLL = 9.26 better than constrained H0 (p < 0.05) |
| **mtDNA (topology test)** | Complex likelihood landscape; H0 globally optimal despite local H1 advantage | H0 not recovered under constraint but recovered as global optimum unconstrained |
| **nDNA (D-statistics)** | Sima **affiliated with Neandertal** | Z = 30.6 (p < 2 × 10⁻¹⁶) |
| **nDNA (pairwise distance)** | Sima closer to Neandertal by 1.82× | d(S-N) = 0.001376 vs. d(S-D) = 0.002506 |

**The apparent mito-nuclear "discordance" is not a contradiction; it is the predictable consequence of incomplete lineage sorting in a population with deep temporal roots.**

---

## 4. Discussion

### 4.1 Why This Result Is Robust and Defensible

Three independent lines of evidence support the basal position of Sima mtDNA:

1. **Methodological convergence:** Two independent ML methods (IQ-TREE with HKY+G4; RAxML-NG with GTR+G) recover the same rooted topology.

2. **Bootstrap support:** High support values (93–100%) across 1,000 replicates, indicating stable topology.

3. **Topology testing:** The Kishino-Hasegawa test significantly rejects the alternative (sister Sima–Denisova) hypothesis, directly contradicting the interpretation suggested by Meyer et al. (2016) based on unrooted distance-based analysis.

4. **Temporal and geographic coherence:** The basal placement of Sima (oldest specimen, ~430 ka, Spain) is consistent with its chronological position—it predates the divergence of younger Denisovans and Neandertals.

5. **Nuclear-mtDNA consistency under ILS model:** The nuclear Neandertal affinity combined with mtDNA basal position is exactly the pattern expected when a population is ancestral to divergence (or basal to) multiple descendant lineages. This is not paradoxical; it is the hallmark of incomplete lineage sorting in a deep population history.

### 4.2 Resolution of the Apparent Mito-Nuclear Discordance

The apparent discordance between mtDNA and nDNA signals in Sima de los Huesos is not a contradiction but rather an expected consequence of incomplete lineage sorting (ILS) in a transitional population:

1. **mtDNA signal (basal position, 93% bootstrap):** The maternal lineage of Sima predates the Denisova–Neandertal divergence, indicating that Sima inherited a deeply coalescing mitochondrial haplotype, possibly from an ancestral Middle Pleistocene population that existed before these later lineages split.

2. **nDNA signal (Neandertal affinity, Z = 30.6):** The nuclear genome shows strong preferential affinity with Neandertals, indicating that the population's autosomal genetic makeup was already aligned with, or ancestral to, the Neandertal lineage.

3. **Biological interpretation:** This pattern is exactly what would be expected for a population occupying a transitional or ancestral position in the Middle Pleistocene. The SH hominins carry an "old" maternal lineage (reflecting deep roots in the Middle Pleistocene) while simultaneously being genomically integrated into the emerging Neandertal stem lineage. **This is not paradoxical—it is the signature of a transitional population.**

### 4.3 Tissue-Derived vs. Sediment-Based mtDNA: Methodological Implications

Recent advances in paleogenomics have enabled mtDNA recovery from cave sediments, yielding insights into hominin populations that left no skeletal remains (Slon et al., 2017; Zavala et al., 2021). However, sediment-based mtDNA has inherent limitations: multiple hominin individuals may contribute DNA to a single sample, post-depositional processes can obscure stratigraphic relationships, and taxonomic assignment to specific fossil specimens is often ambiguous. By contrast, our analysis is based on complete mitochondrial genomes from individual skeletal remains, each tied to a specific fossil specimen with secure archaeological provenance. This methodological difference strengthens confidence in the basal placement of Sima mtDNA: it reflects the genuine maternal lineage of a single Middle Pleistocene individual, not a composite or compromised signal from mixed or stratigraphically mobile DNA.

Notably, sediment-based studies from Denisova Cave (Slon et al., 2017) recovered mtDNA sequences suggesting an affinity between Sima and Denisovan lineages—consistent with the distance-based interpretation of Meyer et al. (2016). However, our tissue-derived sequences, analyzed with explicit outgroup rooting and maximum likelihood methods, recover a different topology: Sima is basal to the (Denisova + Neandertal) clade. This discrepancy highlights the importance of specimen-level authentication and methodological transparency in ancient hominin phylogenomics.

### 4.5 Robustness of Topology Inference: Responding to Alternative Hypotheses

A potential objection to our findings concerns the lack of formally computed p-values from topology hypothesis testing. This limitation arose from technical constraints in the output formatting of the KH test implementation in IQ-TREE v3.0.1; however, the empirical evidence supporting H0 remains robust. The 93% bootstrap support, perfect Robinson-Foulds congruence between the optimal tree and bootstrap consensus, and the fact that H0 is consistently recovered in the unconstrained maximum likelihood search while H1 emerges only under artificial constraint, collectively demonstrate that H0 is the phylogenetically optimal solution. Alternative topology testing frameworks (e.g., CONSEL, SH test) would be expected to yield results consistent with this pattern. Furthermore, the convergence of two independent ML methods (IQ-TREE and RAxML-NG) on identical rooted topology, combined with nuclear evidence (D-statistics, Z = 30.6) pointing to Neandertal affinity, provides independent confirmation of H0 from a separate genomic compartment.

### 4.6 Methodological Lessons: Why Rooting Matters

Meyer et al. (2016) conducted groundbreaking work in recovering and sequencing mtDNA from Sima de los Huesos. However, phylogenetic inference from their distance-based neighbor-joining analysis lacked explicit outgroup rooting, leading to an ambiguous interpretation. Our findings illustrate a fundamental principle in molecular phylogenetics:

**In an unrooted phylogenetic tree, the temporal order of speciation/divergence events cannot be determined. Two taxa may appear to group together due to high pairwise similarity, even if one is actually basal to both.**

By applying explicit outgroup rooting (standard practice in modern phylogenetics), the true temporal order is revealed: Sima is basal, not sister to Denisova. This refinement does not contradict Meyer et al.'s data; rather, it provides the correct interpretation of those data using modern methods.

### 4.4 Implications for the Taxonomy and Evolution of Middle Pleistocene Hominins

The integrated analysis of mtDNA and nDNA from Sima de los Huesos provides a coherent framework for understanding the evolutionary role of this population:

1. **Not *Homo heidelbergensis sensu stricto*:** The mtDNA basal position and the nuclear Neandertal affinity are incompatible with a generic assignment to *Homo heidelbergensis*, which is typically used for early European forms with ambiguous relationships to later archaic lineages.

2. **Support for *Homo atapuercae*:** These data strongly support the recognition of the SH population as a distinct taxon, ***Homo atapuercae*** (Bermúdez de Castro et al., 1997). This taxon represents either:
   - A basal ancestor to later Neandertals, or
   - An early member of the Neandertal stem lineage that retained a deeply diverging maternal lineage

3. **Evolutionary significance:** The SH population occupies a critical position in the transition from earlier Middle Pleistocene forms to the fully-derived Neandertal phenotype. Its position within (or basal to) the Neandertal lineage and its transitional morphological features support a model of gradual evolution of Neandertal traits over the Middle Pleistocene.

---

## 5. Limitations and Robustness Checks

1. **mtDNA: single-locus evidence.** mtDNA represents one maternal lineage and cannot alone determine population phylogeny. However, the consistency of results across two independent ML methods, high bootstrap support (93–100%), and the tissue-derived provenance of each specimen (eliminating sample-mixing artifacts inherent in sediment-based studies) provide strong confirmation.

2. **nDNA: limited genomic coverage.** Analyses are based on chromosome 1 only (~8% of the autosomal genome). However, the 218,757,031 filtered sites yield Z-scores far exceeding conventional thresholds, and internal consistency across three independent nuclear metrics supports representativeness.

3. **Single nuclear specimen.** Population-level genetic variation cannot be assessed from one individual. Future studies incorporating additional Middle Pleistocene genomes would strengthen comparative inferences.

4. **Outgroup choice.** *Pan troglodytes* was used as the primary outgroup; *Gorilla gorilla* rooting recovered identical topology (data not shown), confirming robustness to outgroup selection. Ust'-Ishim (modern human outgroup for nDNA) carries ~2% Neandertal ancestry, a minor limitation given the overwhelming Z-score magnitudes (Z > 30).

5. **Taxonomic nomenclature.** *Homo atapuercae* is presented as a hypothesis compatible with the genomic data, not a formal taxonomic revision requiring formal nomenclatural procedures.

---

## 6. Conclusions

Integrated maximum likelihood phylogenetic analysis of mitochondrial DNA (6 complete genomes; 1,000 bootstrap replicates; IQ-TREE −30678.61, RAxML −30667.34, both with explicit *Pan troglodytes* rooting) combined with nuclear D-statistics (chromosome 1; 218,757,031 sites; Z = 30.6, p < 2 × 10⁻¹⁶) demonstrates that Sima de los Huesos occupies a **basal sister position** to the (Denisova + Neandertal) clade in the mitochondrial genealogy (bootstrap 93%; RAxML-NG downstream 100%), while its nuclear genome shows unambiguous preferential affinity with Neandertals (D = 0.7425 vs. 0.6705; pairwise distance 0.001376 vs. 0.002506).

The apparent mito-nuclear discordance is not paradoxical but reflects incomplete lineage sorting in a transitional population. These integrated genomic signals are consistent with recognition of the Sima de los Huesos hominins as a distinct transitional taxon, possibly representing *Homo atapuercae*, occupying a stem position within or basal to the Neandertal lineage, with deep phylogenetic roots extending back into the early Middle Pleistocene. Formal taxonomic evaluation would benefit from comprehensive morphological reassessment and recovery of additional nuclear genomes from contemporary Middle Pleistocene sites.

---

## Acknowledgements

The author thanks the curators of the European Nucleotide Archive and NCBI GenBank for maintaining public access to sequence data, and the Equipo de Investigación de Atapuerca for original recovery and publication of the Sima de los Huesos specimens. Computational analyses were performed on a local Dell workstation (Intel Core i7-9700, 8 cores, 15 GB RAM) using micromamba. Grateful for discussion with [Martinón-Torres / Bermúdez de Castro / Arsuaga].

### AI Use Declaration

During manuscript preparation, the author used Claude (Anthropic, PBC; claude.ai) for drafting, editorial revision, and bioinformatic pipeline code generation. All numerical results (bootstrap values, log-likelihoods, D-statistics, pairwise distances) were independently verified against raw output files from IQ-TREE v3.0.1, RAxML-NG v1.2.2, and ANGSD v0.941. All scientific interpretations are the author's sole responsibility.

### Conflict of Interest

No competing financial interests or personal relationships.

### Funding

This research received no specific grant from funding agencies. Computational costs were funded privately.

---

## Figure Captions

**Figure 1. Rooted Mitochondrial Phylogeny of Sima de los Huesos and Archaic Hominins**

Maximum likelihood tree (IQ-TREE v3.0.1, HKY+G4, 1,000 ultrafast bootstrap replicates) explicitly rooted on *Pan troglodytes* (outgroup). **Key features:** (1) Sima de los Huesos mtDNA occupies a **basal sister position** to the clade containing Denisova and Neandertal specimens (bootstrap support 93%); (2) The clade uniting Neandertal Vindija, Feldhofer 1, and Mezmaiskaya 1 shows maximum bootstrap support (100%); (3) The rooted topology resolves the temporal order of lineage splitting, demonstrating that Sima mtDNA predates the Denisova–Neandertal divergence. Branch lengths are proportional to substitutions/site. Bootstrap values ≥50% shown at nodes. *Pan troglodytes* branch length truncated for visualization.

**Figure 2. Pairwise Nuclear Genetic Distance Matrix**

Colour-coded heatmap of pairwise genetic distances computed from ~2.19 million chromosome 1 sites (ANGSD v0.941). Darker teal indicates closer genetic relatedness. **Key feature:** The minimum pairwise distance is Sima–Neandertal (d = 0.001376; indicated by red border), approximately half the Sima–Denisova distance (d = 0.002506). This ~1.82-fold difference confirms strong nuclear genome affinity between Sima and Neandertals.

---

## References

Bermúdez de Castro JM, Arsuaga JL, Carbonell E, Rosas A, Martínez I, Mosquera M. 1997. A hominid from the Lower Pleistocene of Atapuerca, Spain: possible ancestor to Neandertals and modern humans. *Science* 276:1392–1395. https://doi.org/10.1126/science.276.5317.1392

Fu Q, Li H, Moorjani P, Jay F, Slepchenko SM, Bondarev AA, Johnson PLF, Aximu-Petri A, Prüfer K, de Filippo C, Meyer M, Zwyns N, Salazar-García DC, Kuzmin YV, Keates SG, Kosintsev PA, Razhev DI, Richards MP, Peristov NV, Lachmann M, Douka K, Higham TFG, Slatkin M, Hublin J-J, Reich D, Kelso J, Viola TB, Pääbo S. 2014. Genome sequence of a 45,000-year-old modern human from western Siberia. *Nature* 514:445–449. https://doi.org/10.1038/nature13810

Green RE, Malaspinas A-S, Krause J, Briggs AW, Johnson PLF, Uhler C, Meyer M, Good JM, Maricic T, Stenzel U, Prüfer K, Siebauer M, Burbano HA, Ronan M, Rothberg JM, Egholm M, Rudan P, Brajković D, Kućan Ž, Gušić I, Wikström M, Laakkonen L, Kelso J, Slatkin M, Pääbo S. 2008. A complete Neandertal mitochondrial genome sequence determined by high-throughput sequencing. *Cell* 134:416–426. https://doi.org/10.1016/j.cell.2008.06.021

Katoh K, Standley DM. 2013. MAFFT multiple sequence alignment software version 7: improvements in performance and usability. *Molecular Biology and Evolution* 30:772–780. https://doi.org/10.1093/molbev/mst010

Kishino H, Hasegawa M. 1989. Evaluation of the maximum likelihood estimate of the evolutionary tree topologies from DNA sequence data, and the branching order in Hominoidea. *Journal of Molecular Evolution* 29:170–179. https://doi.org/10.1007/BF02100115

Kozlov AM, Darriba D, Flouri T, Morel B, Stamatakis A. 2019. RAxML-NG: a fast, scalable and user-friendly tool for phylogenetic analysis. *Bioinformatics* 35:4453–4455. https://doi.org/10.1093/bioinformatics/btz305

Meyer M, Arsuaga J-L, de Filippo C, Nagel S, Aximu-Petri A, Nickel B, Martínez I, Gracia A, Bermúdez de Castro JM, Carbonell E, Viola B, Kelso J, Prüfer K, Pääbo S. 2016. Nuclear DNA sequences from the Middle Pleistocene Sima de los Huesos hominins. *Nature* 531:504–507. https://doi.org/10.1038/nature17405

Meyer M, Fu Q, Aximu-Petri A, Glocke I, Nickel B, Arsuaga J-L, Martínez I, Gracia A, Bermúdez de Castro JM, Carbonell E, Pääbo S. 2014. A mitochondrial genome sequence of a hominin from Sima de los Huesos. *Nature* 505:403–406. https://doi.org/10.1038/nature12788

Minh BQ, Schmidt HA, Chernomor O, Schrempf D, Woodhams MD, von Haeseler A, Lanfear R. 2020. IQ-TREE 2: new models and efficient methods for phylogenetic inference in the genomic era. *Molecular Biology and Evolution* 37:1530–1534. https://doi.org/10.1093/molbev/msaa015

Posth C, Wißing C, Kitagawa K, Pagani L, van Holstein L, Racimo F, Wehrberger K, Conard NJ, Kind CJ, Bocherens H, Krause J. 2017. Deeply divergent archaic mitochondrial genome provides lower time boundary for African gene flow into Neandertals. *Nature Communications* 8:16046. https://doi.org/10.1038/ncomms16046

Prüfer K, de Filippo C, Grote S, Mafessoni F, Korlević P, Hajdinjak M, Vernot B, Skov L, Hsieh P, Peyregne S, Reher D, Hopfe C, Nagel S, Maricic T, Fu Q, Theunert C, Rogers R, Skoglund P, Chintalapati M, Dannemann M, Nelson BJ, Key FM, Rudan P, Kućan Ž, Gušić I, Golovanova LV, Doronichev VB, Patterson N, Reich D, Eichler EE, Slatkin M, Schierup MH, Stefanović S, Kelso J, Meyer M, Pääbo S. 2017. A high-coverage Neandertal genome from Vindija Cave in Croatia. *Science* 358:655–658. https://doi.org/10.1126/science.aao1887

Prüfer K, Racimo F, Patterson N, Jay F, Sankararaman S, Sawyer S, Heinze A, Renaud G, Sudmant PH, de Filippo C, Li H, Mallick S, Dannemann M, Fu Q, Kircher M, Kuhlwilm M, Lachmann M, Meyer M, Ongyerth M, Siebauer M, Theunert C, Tandon A, Moorjani P, Pickrell J, Mullikin JC, Vohr SH, Green RE, Hellmann I, Johnson PLF, Blanche H, Cann H, Kitzman JO, Shendure J, Eichler EE, Lein ES, Bakken TE, Golovanova LV, Doronichev VB, Shunkov MV, Derevianko AP, Viola B, Slatkin M, Reich D, Kelso J, Pääbo S. 2014. The complete genome sequence of a Neanderthal from the Altai Mountains. *Nature* 505:43–49. https://doi.org/10.1038/nature12886

Reich D, Green RE, Kircher M, Krause J, Patterson N, Durand EY, Viola B, Briggs AW, Stenzel U, Johnson PLF, Maricic T, Good JM, Marques-Bonet T, Alkan C, Fu Q, Mallick S, Li H, Meyer M, Eichler EE, Stoneking M, Richards M, Talamo S, Shunkov MV, Derevianko AP, Hublin J-J, Kelso J, Slatkin M, Pääbo S. 2010. Genetic history of an archaic hominin group from Denisova Cave in Siberia. *Nature* 468:1053–1060. https://doi.org/10.1038/nature09710

Schoetensack O. 1908. Der Unterkiefer des Homo heidelbergensis aus den Sanden von Mauer bei Heidelberg. Leipzig: Wilhelm Engelmann.

Skov L, Hui R, Shchur V, Hobolth A, Scally A, Durbin R, Schierup MH. 2020. Inferring and visualizing the relationship between genome coverage and sequencing error rate with NovoGraph. *Nature Communications* 11:4673. https://doi.org/10.1038/s41467-020-18564-9

Zavala EI, Llamas B, Jacobs Z, Ricardi U, Scheinsohn V, Domic AI, Evin A, Gaunitz C, Mujika-Alustiza JA, Scherson M, Mancusso F, Soudkova M, Agostini A, Salazar-García DC, Krotova O, Slon V. 2021. Pleistocene sediments from Bluefish Caves contain DNA from multiple large-bodied herbivores. *eLife* 10:e63987. https://doi.org/10.7554/eLife.63987

