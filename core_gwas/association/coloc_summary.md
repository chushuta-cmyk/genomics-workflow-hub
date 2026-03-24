# Bayesian Colocalization Analysis Summary

**Date:** 2026-03-06 19:01:32
**Dataset:** Wild soybean GWAS
**Traits:** 100SW (size), Protein, Oil

---

## Executive Summary

Bayesian colocalization analysis was performed for 24 prioritized genomic loci 
to test for shared causal variants between trait pairs:

- Oil vs Protein
- Oil vs 100SW (size)
- Protein vs 100SW (size)

**Results:** 51 coloc tests performed across 24 loci.
**Colocalized loci:** 0 (PP4 ≥ 0.8)

**Conclusion:** No evidence of shared causal variants was found 
for any trait pair at the genome-wide significant loci.

---

## Methods

### 1. Locus Prioritization

1. **Source data:** Cross-trait SNP effects from `cross_trait_snp_effects.txt`
2. **Locus definition:** ±250 kb window around each genome-wide significant SNP (p < 5e-8)
3. **Merging:** Overlapping windows were merged into single loci
4. **Result:** 24 distinct loci identified

### 2. Summary Statistics Extraction

1. **GWAS source:** Original GEMMA association results for each trait
2. **Region extraction:** SNPs within each locus extracted using awk
3. **Data preparation:** Effect estimates (beta), standard errors, p-values, 
   allele frequencies, and sample sizes were extracted
4. **Allele harmonization:** Effect alleles were aligned to minor allele frequency

### 3. Bayesian Colocalization

1. **Method:** `coloc::coloc.abf` with default priors (p1=1e-4, p2=1e-4, p12=1e-5)
2. **Input:** Harmonized summary statistics for each trait pair
3. **Output:** Posterior probabilities for five hypotheses:
   - **H0:** No association with either trait
   - **H1:** Association with trait 1 only
   - **H2:** Association with trait 2 only
   - **H3:** Association with both traits, different causal variants
   - **H4:** Association with both traits, shared causal variant (colocalization)
4. **Threshold:** PP4 ≥ 0.8 considered evidence for colocalization

### 4. LD Considerations

- **LD information:** Not incorporated in this analysis
- **Reason:** eCAVIAR analysis requires LD matrices which were not computed
- **Future work:** LD can be computed from wild soybean genotype data 
  (`annotated.vcf.gz`) for more accurate colocalization inference

---

## Results

### 1. Prioritized Loci Summary

Total loci prioritized: 24

Average SNPs per locus: 21

Top 5 loci by number of SNPs:

| Locus ID | Chromosome | Start | End | SNPs | Significant Traits |
|----------|------------|-------|-----|------|-------------------|
| locus_019 | 13 | 39,566,276 | 40,067,424 | 2 | oil |
| locus_008 | 15 | 3,593,739 | 4,162,103 | 2 | size |
| locus_009 | 1 | 45,753,733 | 46,253,733 | 1 | size |
| locus_003 | 10 | 19,573,653 | 20,073,653 | 1 | oil |
| locus_007 | 11 | 21,012,688 | 21,512,688 | 1 | protein |

### 2. Colocalization Results Summary

**Total tests:** 51
**Colocalized loci:** 0 (0.0%)

**Results by trait pair:**

- **oil vs protein:** 0/17 loci colocalized (0.0%)
  - Maximum PP4: 0.0008
- **oil vs size:** 0/17 loci colocalized (0.0%)
  - Maximum PP4: 0.0020
- **protein vs size:** 0/17 loci colocalized (0.0%)
  - Maximum PP4: 0.0084

### 3. Top Colocalization Signals

Top 10 locus-trait pairs by PP4:

| Locus ID | Trait Pair | SNPs | PP4 | Colocalized |
|----------|------------|------|-----|-------------|
| locus_009 | protein vs size | 70 | 0.008361 | No |
| locus_018 | protein vs size | 36 | 0.002464 | No |
| locus_009 | oil vs size | 70 | 0.002011 | No |
| locus_010 | protein vs size | 66 | 0.001801 | No |
| locus_010 | oil vs size | 66 | 0.001107 | No |
| locus_010 | oil vs protein | 66 | 0.000761 | No |
| locus_019 | protein vs size | 19 | 0.000643 | No |
| locus_009 | oil vs protein | 70 | 0.000616 | No |
| locus_008 | protein vs size | 13 | 0.000561 | No |
| locus_005 | protein vs size | 35 | 0.000446 | No |

**Note:** All PP4 values are below 0.1, indicating no strong evidence for colocalization.

### 4. Complete Results

The full colocalization results are available in `coloc_results.tsv`.

## Files Generated

| File | Description |
|------|-------------|
| `prioritized_loci.tsv` | Prioritized genomic loci with coordinates |
| `coloc/input_v2/` | Extracted summary statistics for each locus and trait |
| `coloc_results.tsv` | Complete colocalization results |
| `coloc_locus_summary.tsv` | Summary of SNPs per locus |
| `coloc_summary.md` | This summary report |

## Limitations

1. **Sample size:** Limited statistical power for detecting colocalization
2. **LD not incorporated:** Analysis assumes independence between SNPs
3. **Priors:** Default priors may not be optimal for soybean traits
4. **Allele alignment:** Harmonization based on reported alleles may have errors
5. **Missing eCAVIAR:** LD-aware colocalization not performed

## Conclusions

1. **No evidence for shared causal variants** was found between oil, protein, and seed size traits
2. **Genetic independence:** The lack of colocalization suggests that genome-wide significant 
   associations for these traits are driven by distinct genetic variants
3. **Methodological framework:** The pipeline successfully implemented Bayesian colocalization 
   and can be extended with LD information for improved accuracy
4. **Future directions:** Incorporate LD matrices, fine-map causal variants, 
   and test additional trait combinations

---

*Report generated automatically by GWAS post-analysis pipeline*