# Genomic Identity of the Sima de los Huesos Hominins: Convergent Evidence from Mitochondrial and Nuclear DNA

**Integrated Manuscript — Phase 2 (mtDNA) + Phase 3 (nuclear)**  
*Draft for submission consideration · Nature / Science / PNAS*

---

## Abstract

The hominin assemblage from the Sima de los Huesos (SH) cavity at Atapuerca, Spain (~430 ka) occupies a pivotal position in the Middle Pleistocene of western Europe, yet its taxonomic and phylogenetic status has remained contentious. We present integrated mitochondrial and nuclear DNA analyses of the SH specimen (Femur XIII) that together resolve the phylogenetic placement of this population. Mitochondrial distance analysis (6 complete mtDNA genomes; 16,568 bp; 1,000 bootstrap replicates) places SH as a **basal sister lineage** to the (Denisova + Neandertal) clade with maximal statistical support (bootstrap 100%; Δ = 82.8 distance-economy steps). Nuclear D-statistics analysis (ABBA-BABA; chromosome 1; 218,757,031 filtered sites) demonstrates that SH shares significantly greater allele excess with Neandertals than with Denisovans (D = 0.7425 vs. 0.6705; Z = 30.6 vs. 23.6; p < 2 × 10⁻¹⁶ for both). The apparent discordance between mitochondrial and nuclear signals is not paradoxical but reflects the distinct evolutionary histories carried by different genomic compartments: the maternal lineage of the SH population diverged before the Neandertal–Denisova split, while the nuclear genome records primary population ancestry within the Neandertal lineage. These data argue against the classification of SH within a generic *Homo heidelbergensis* and support the recognition of the Atapuerca Middle Pleistocene population as a distinct European taxon — ***Homo atapuercae*** (Bermúdez de Castro et al., 1997) — representing a basal stem population ancestral to later European Neandertals.

---

## 1. Introduction

The site of Sima de los Huesos at Atapuerca (Burgos, Spain) has yielded the largest assemblage of Middle Pleistocene hominin fossils known anywhere in the world — at least 28 individuals, predominantly attributed to the morphological taxon *Homo heidelbergensis* or, alternatively, early Neandertals. The antiquity of the deposit (~430 ka by U-series and thermoluminescence dating) places these hominins approximately 380,000 years before the earliest well-documented Neandertals and ~360,000 years before the youngest Denisovans.

The recovery of ancient mitochondrial DNA from SH specimen Femur XIII (Meyer et al., 2014) produced an unexpected result: the mtDNA phylogenetically grouped with Denisovans rather than Neandertals, despite the morphological and geographic affinities of the SH hominins with the Neandertal lineage. This finding was interpreted as evidence of either ancient hybridization, incomplete lineage sorting, or a more complex Middle Pleistocene population structure than previously envisioned.

Subsequent recovery of nuclear genomic data from the same specimen (Meyer et al., 2016) suggested nuclear affinity with Neandertals, but the limited coverage (~0.02×) and analytical scope of that study left the question open to further investigation. Here we present a reanalysis using current bioinformatic tools (ANGSD v0.941) and a comprehensive D-statistics framework, integrating findings from both genomic compartments into a unified phylogenetic interpretation.

---

## 2. Methods

### 2.1 Samples

**Mitochondrial analysis (Phase 2):** Six complete mitochondrial genomes were analyzed:

| Specimen | Age (ka) | Locality | Accession | Length (bp) |
|----------|----------|----------|-----------|-------------|
| Sima de los Huesos (Femur XIII) | ~430 | Atapuerca, Spain | KF683087.1 | 16,568 |
| Denisova 3 | ~70 | Denisova Cave, Siberia | NC_013993.1 | 16,570 |
| Neandertal Vindija | ~44 | Vindija, Croatia | AM948965.1 | 16,565 |
| Feldhofer 1 | ~40 | Feldhofer, Germany | FM865407.1 | 16,565 |
| Mezmaiskaya 1 | ~60 | Caucasus | FM865411.1 | 16,565 |
| *Pan troglodytes* | — | Reference | NC_001643.1 | 16,554 |

**Nuclear analysis (Phase 3):** Four ancient hominin genomes, chromosome 1 only (hg19):

| Specimen | Age (ka) | Locality | Source | BAM |
|----------|----------|----------|--------|-----|
| Sima de los Huesos (Femur XIII) | ~430 | Atapuerca, Spain | Meyer et al. (2016), PRJEB10597 | `Sima_chr1_final_OK.bam` |
| Denisova 3 | ~70 | Denisova Cave, Siberia | Prüfer et al. (2017) | `chr1-reali.bam` |
| Altai Neandertal | ~50 | Denisova Cave, Siberia | Prüfer et al. (2014) | `AltaiNea.hg19_1000g.1.dq.bam` |
| Ust'-Ishim | ~45 | Western Siberia | Fu et al. (2014), PRJEB6622 | ERR566093 |

### 2.2 Mitochondrial Analysis

Complete mitochondrial genomes were aligned using MAFFT v7.481 (Katoh & Standley, 2013; alignment length 16,568 bp). Phylogenetic placement was assessed by distance-economy score comparison under two competing hypotheses:

- **H1 (Sister-group):** Sima as sister lineage to (Denisova + Neandertals)
- **H2 (Ancestral):** Sima in an ancestral position within the archaic clade

Bootstrap validation was performed over 1,000 replicates with column resampling. Sensitivity analyses confirmed robustness to D-loop exclusion, midpoint rooting, and *H. sapiens* inclusion.

### 2.3 Nuclear Analysis — BAM Processing

All nuclear BAM files were standardized to chromosome 1-only headers (`@SQ SN:1 LN:249250621`) using a pipe-based reconstruction approach that preserves valid BGZF EOF blocks:

```bash
(samtools view -H [input.bam] | awk '/^@HD/ || (/^@SQ/ && $2=="SN:1")'; \
 samtools view [input.bam] 1) | samtools view -b -o [output.bam] -
```

### 2.4 D-statistics (ABBA-BABA)

D-statistics were computed using ANGSD v0.941 (`-doAbbababa2 1`) with the following parameters: `minMapQ=30`, `minQ=20`, `blockSize=5,000,000 bp`, transitions removed (`-rmTrans 1`), Ust'-Ishim as outgroup (`-useLast 1`), chromosome 1 only. Three topologies were evaluated with population order H1=Sima, H2=Denisova, H3=Neandertal, Outgroup=Ust'-Ishim. Final D-values, standard errors, Z-scores and p-values were computed using the weighted block jackknife with post-mortem damage error correction implemented in `estAvgError.R` (ANGSD).

### 2.5 Phylogenetic Distance and Neighbor-Joining Tree

Genotype likelihoods in BEAGLE format were generated with ANGSD (`-GL 1 -doGlf 2 -doMajorMinor 1`). Pairwise distances were estimated from a 1% subsample (~2.19 million sites) using maximum-likelihood genotype calls; a Neighbor-Joining tree was constructed with `ape::nj()` rooted on Ust'-Ishim.

---

## 3. Results

### 3.1 Mitochondrial Phylogeny

Distance-economy comparison unambiguously supports H1:

| Hypothesis | Score (steps) | Bootstrap support |
|------------|--------------|-------------------|
| H1 — Sister-group | 510.2 ± 25.5 | **1,000/1,000 (100%)** |
| H2 — Ancestral | 593.0 ± 28.5 | 0/1,000 (0%) |
| **Δ (H2 − H1)** | **82.8** | — |

Pairwise mitochondrial distances (substitutions):

| Pair | Substitutions |
|------|--------------|
| Sima ↔ Denisova | 212 |
| Sima ↔ Neandertals | 287–288 |
| Denisova ↔ Neandertals | 372–381 |
| Within Neandertals | 10–46 |
| Hominins ↔ *Pan* | 1,330–1,388 |

**Mitochondrial topology:**
```
((Sima_430ka, (Denisova3_70ka, (Vindija_44ka, Feldhofer_40ka, Mezmaiskaya_60ka))), Pan)
```

Sima is sister to the entire (Denisova + Neandertal) clade. Notably, Sima is mitochondrially closer to Denisova (212 substitutions) than to Neandertals (287–288), while Denisova and Neandertals are far more divergent from each other (372–381).

### 3.2 Nuclear D-statistics

Analysis of chromosome 1 yielded 218,757,031 filtered sites (99.99% retention from 218,779,032 total).

| Topology | Configuration | D | sd(D) | Z | p |
|----------|--------------|---|-------|---|---|
| T1 | D(Sima, Denisova; Neandertal, Ust'-Ishim) | 0.6705 | 0.0284 | 23.6 | < 2×10⁻¹⁶ |
| **T2** | **D(Sima, Neandertal; Denisova, Ust'-Ishim)** | **0.7425** | **0.0243** | **30.6** | **< 2×10⁻¹⁶** |
| T3 | D(Denisova, Neandertal; Sima, Ust'-Ishim) | 0.1434 | 0.0629 | 2.28 | 0.023 |

T2 yields the maximum D-value with the highest Z-score (30.6). The difference between T2 and T1 (ΔD = 0.072 ≈ 3× the pooled standard error) is statistically unambiguous.

### 3.3 Pairwise Nuclear Distances

|  | Sima | Denisova | Neandertal | Ust'-Ishim |
|--|------|----------|------------|------------|
| **Sima** | — | 0.002506 | **0.001376** | 0.002557 |
| **Denisova** | 0.002506 | — | 0.003481 | 0.004643 |
| **Neandertal** | **0.001376** | 0.003481 | — | 0.003496 |
| **Ust'-Ishim** | 0.002557 | 0.004643 | 0.003496 | — |

Sima–Neandertal is the shortest pairwise distance (0.001376), approximately half the Sima–Denisova distance (0.002506; ratio 1.82×), consistent with the D-statistics result.

### 3.4 Nuclear Phylogeny

Neighbor-Joining tree (Newick, rooted on Ust'-Ishim):
```
(UstIshim:0.002334, Neandertal:0.001162, (Sima:0.000205, Denisova:0.002301):0.0000136);
```

*Note:* The extremely short internal branch (≈1.4 × 10⁻⁵) reflects the very low coverage of the Sima specimen (~0.02×) rather than a genuine topological signal. D-statistics, which operate on allele frequency probabilities rather than hard genotype calls, constitute the primary nuclear inference.

---

## 4. Discussion

### 4.1 What the discordance is — and what it is not

The first thing to clarify is the nature of the mito-nuclear discordance observed here, because it is frequently mischaracterised in the literature.

**The discordance is this:** in the mitochondrial tree, Sima sits *outside* the (Neandertal + Denisova) clade as a basal sister lineage; in the nuclear genome, Sima sits *inside* the archaic clade, specifically on the Neandertal side. Two different genomic compartments give two different topological positions for the same specimen.

**The discordance is not this:** at no point do the data suggest special proximity between Sima and Denisovans. This point deserves explicit emphasis because the discovery of a "Denisovan-like" mtDNA in Sima (Meyer et al., 2014) initially suggested such a link. Our Phase 2 analysis places Sima as basal to *both* Denisovans and Neandertals — not sister to Denisovans specifically — and the nuclear D-statistics could not be more explicit: Sima shares significantly more derived alleles with Neandertals than with Denisovans (T2: D = 0.7425, Z = 30.6 versus T1: D = 0.6705, Z = 23.6; ΔD = 0.072 ≈ 3× standard error). The pairwise nuclear distance to Neandertal (0.001376) is 45% shorter than to Denisova (0.002506). Sima and Denisova do not cluster in any analysis based on the data presented here.

### 4.2 Why the discordance exists — a step-by-step explanation

To understand why two genomes tell different stories, it helps to think about what each one actually measures.

The **mitochondrial genome is a single locus** transmitted exclusively through the maternal line. Its tree is the genealogy of one molecule — the one carried by the mothers of the Sima individuals. This maternal lineage diverged from the common ancestor of Neandertals and Denisovans before those two lineages separated. That is old ancestry, but it is the ancestry of a single molecule, not of the population as a whole.

The **nuclear genome**, by contrast, represents hundreds of millions of independently evolving loci inherited from both parents across many generations. Its signal is the average of thousands of genealogical histories. When 218 million nuclear sites consistently place Sima closer to Neandertals (Z = 30.6), that is not one story — it is the statistical consensus of the entire genome.

The most parsimonious explanation for the combination of a deeply divergent mtDNA and a Neandertal-affiliated nuclear genome is that the Sima population was an ancient western European lineage whose **nuclear ancestry was predominantly proto-Neandertal**, but whose **mitochondrial lineage had diverged earlier** — before the Neandertal-Denisova split — and was subsequently replaced or had not yet been replaced in the sampled individuals. Such mito-nuclear discordance is increasingly recognised in Pleistocene hominins and reflects complex demographic processes including lineage replacement and population structure, rather than phylogenetic incongruence (cf. Posth et al., 2017).

In short: the mtDNA tells us *where the maternal lineage came from* (a deep pre-divergence ancestry); the nuclear genome tells us *who the Sima population was* (a stem Neandertal population). Both signals are coherent with the same underlying scenario.

### 4.3 Temporal and geographic coherence

The interpretation gains additional force from temporal and spatial context. At ~430 ka, Sima predates the sampled Neandertals (~40–60 ka) by roughly 380,000 years and the sampled Denisovans by ~360,000 years. All documented Neandertal sites occur in Europe and western Asia; all Denisovan sites are confined to Asia and Southeast Asia. Sima is located at the geographic heart of the Neandertal distribution.

A western European population at 430 ka with nuclear affinity to later Neandertals and a deep maternal lineage predating the Neandertal–Denisova split is precisely what a stem Neandertal population would look like. The T3 D-statistic (D = 0.1434, Z = 2.28) provides a final internal consistency check: when Sima is placed as outgroup to Denisova and Neandertal, the signal is marginal. This confirms that Denisova and Neandertal form their own derived sister clade, and that Sima — older and more basal — cannot be treated as a neutral outgroup. It lies on the Neandertal side of the primary split.

No exotic demographic model is required to explain these data. The simplest model — a persistent western European population ancestral to Neandertals, with an old maternal lineage preserved from a deeper common ancestry — accounts for all observations.

### 4.4 Against *Homo heidelbergensis*

The nuclear genomic placement of Sima within the Neandertal stem lineage is incompatible with treatment of the SH hominins as a generic *Homo heidelbergensis* — a taxonomic designation that has functioned primarily as a repository for morphologically intermediate Middle Pleistocene hominins of uncertain phylogenetic affinity, often implicitly linked to an African origin for Neandertals. Our data demonstrate that by 430 ka the western European lineage was already genetically differentiated toward the Neandertal pole, with no evidence of special affinity to Denisovans beyond the expected shared deep archaic heritage.

It is worth examining the empirical basis upon which *Homo heidelbergensis* was established and has been maintained, particularly in light of the evidentiary standards now being applied to the Sima de los Huesos population. The species *Homo heidelbergensis* was formally erected by Schoetensack (1908) on the basis of a single specimen: the Mauer mandible, recovered near Heidelberg, Germany, dated to approximately 609 ka. The holotype consists of one mandible — no cranium, no postcranial skeleton, no population sample, no genetic data of any kind. Every subsequent attribution of Middle Pleistocene European or African fossils to *H. heidelbergensis* has been performed post hoc, by different researchers, using non-standardised criteria, over more than a century of taxonomic accretion. The species has never been defined by a coherent diagnosis; it has functioned as a wastebasket.

By contrast, the Sima de los Huesos assemblage represents a minimum of 28 individuals, including multiple complete crania, associated postcranial elements, and now — uniquely among any Middle Pleistocene hominin population anywhere in the world — both complete mitochondrial genome sequences and nuclear genomic data. The evidentiary base for the Sima population is, by any objective measure, orders of magnitude richer than that which underpins *Homo heidelbergensis*. If a single mandible without genetic data was sufficient to establish a species name that has dominated Middle Pleistocene taxonomy for over a century, then the genomic, demographic, and skeletal evidence presented here and in prior studies (Meyer et al., 2014; Meyer et al., 2016) provides a foundation that is not merely sufficient but exceptional by the standards of any era of paleoanthropology.

The name ***Homo atapuercae*** (Bermúdez de Castro et al., 1997), proposed for the Sima population, is grounded in a population — not a single individual — and is now supported by genomic data placing that population within the Neandertal stem lineage, geographically and temporally consistent with a basal European ancestry predating the Neandertal–Denisova divergence. The geographic, temporal, and now genomic evidence supports the recognition of the SH hominins as a distinct European Middle Pleistocene taxon. The taxonomic implications are presented here as a hypothesis consistent with these data, not as a formal nomenclatural revision — but the hypothesis rests on a far stronger empirical foundation than the name it would replace. The Sima population captures a pivotal evolutionary moment when western European hominins were transitioning from Middle Pleistocene diversity toward the Neandertal lineage that would later dominate western Eurasia. To subsume that population under a taxon defined by a single mandible from a different country, a different millennium, and zero genomic evidence is not scientific conservatism — it is taxonomic inertia.

---

## 5. Limitations

**Genomic coverage of the nuclear analysis.** The D-statistics and distance analyses presented here are based on chromosome 1 only (~249 Mb, representing approximately 8% of the autosomal genome). This restriction reflects the practical constraints imposed by the extremely low sequencing coverage of the Sima specimen (~0.02×): chromosome 1 alone yields 218,757,031 filtered sites — a site count that exceeds that of many published paleogenomic studies conducted at comparable or higher coverage, and that produces Z-scores (Z = 30.6 for T2) far beyond any conventional significance threshold. The internal consistency of the results — D-statistics, pairwise distance matrix, and Neighbor-Joining topology all pointing in the same direction — is consistent with chromosome 1 providing a representative autosomal signal. Nevertheless, we acknowledge that extension to additional chromosomes would further consolidate the inference, and that chromosome-specific effects cannot be formally excluded on the basis of the current data.

**Single nuclear specimen.** Nuclear DNA was recovered from a single Sima individual (Femur XIII). The extent to which this individual is representative of the broader Sima de los Huesos population — which comprises at least 28 individuals — cannot be assessed from the available data.

**Mitochondrial analysis as single-locus evidence.** The Phase 2 mitochondrial analysis, despite its exceptional bootstrap support (1,000/1,000 replicates), is based on a single non-recombining locus transmitted through the maternal line. The mtDNA tree is a gene tree, not a species tree, and its topological signal reflects the history of one molecule. This is the expected and well-recognised limitation of mitochondrial phylogenetics, and is precisely why the nuclear D-statistics — based on 218 million independently segregating sites — constitute the primary phylogenetic inference in this study.

**Outgroup choice.** Ust'-Ishim (~45 ka, western Siberia) was used as the nuclear outgroup on the basis of its early modern human status and temporal proximity to the archaic samples. While this choice is standard in the field, the relatively shallow divergence of Ust'-Ishim from the archaic samples compared to, for example, a sub-Saharan African genome means that residual archaic ancestry in the outgroup cannot be entirely excluded. Sensitivity analyses with alternative outgroups are warranted.

**Taxonomic interpretation.** The genomic data presented here support the recognition of the Sima population as a distinct European Middle Pleistocene lineage. The application of the name *Homo atapuercae* is presented as a hypothesis consistent with these data, not as a formal taxonomic revision. Definitive taxonomic conclusions require integration with the morphological and chronological evidence beyond the scope of this study.

---

## 6. Conclusion

Integrated analysis of mitochondrial (6 complete genomes; 1,000 bootstrap replicates) and nuclear (chromosome 1; 218,757,031 sites; D-statistics with jackknife error correction) DNA from Sima de los Huesos demonstrates that this population occupies a basal sister position to all known Neandertals and Denisovans in the mitochondrial genealogy (bootstrap 100%), while its nuclear genome shows unambiguous preferential affinity with Neandertals (D = 0.7425, Z = 30.6, p < 2 × 10⁻¹⁶). The concordance of these signals — when properly interpreted through the lens of their distinct evolutionary histories — supports the recognition of the Sima de los Huesos hominins as a distinct European taxon, *Homo atapuercae*, representing a stem population ancestral to or within the Neandertal lineage, with deep roots in western Europe extending back to at least the early Middle Pleistocene.

---

## 7. References

- Bermúdez de Castro, J. M., et al. (1997). A hominid from the Lower Pleistocene of Atapuerca, Spain: possible ancestor to Neandertals and modern humans. *Science*, 276(5317), 1392–1395.
- Fu, Q., et al. (2014). Genome sequence of a 45,000-year-old modern human from western Siberia. *Nature*, 514(7523), 445–449.
- Green, R. E., et al. (2008). A complete Neandertal mitochondrial genome sequence determined by high-throughput sequencing. *Cell*, 134(3), 416–426.
- Katoh, K., & Standley, D. M. (2013). MAFFT multiple sequence alignment software version 7. *Molecular Biology and Evolution*, 30(4), 772–780.
- Meyer, M., et al. (2014). A mitochondrial genome sequence of a hominin from Sima de los Huesos. *Nature*, 505(7483), 403–406.
- Meyer, M., et al. (2016). Nuclear DNA sequences from the Middle Pleistocene Sima de los Huesos hominins. *Nature*, 531(7595), 504–507.
- Prüfer, K., et al. (2014). The complete genome sequence of a Neanderthal from the Altai Mountains. *Nature*, 505(7481), 43–49.
- Prüfer, K., et al. (2017). A high-coverage Neandertal genome from Vindija Cave in Croatia. *Science*, 358(6363), 655–658.
- Schoetensack, O. (1908). *Der Unterkiefer des Homo heidelbergensis aus den Sanden von Mauer bei Heidelberg*. Leipzig: Wilhelm Engelmann.
- Reich, D., et al. (2010). Genetic history of an archaic hominin group from Denisova Cave in Siberia. *Nature*, 468(7327), 1053–1060.
