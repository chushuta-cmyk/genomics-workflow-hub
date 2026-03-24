# Wild Soybean GWAS Post-Analysis Final Report

## Executive Summary

This report summarizes the complete post-GWAS analysis of wild soybean traits: 
100-seed weight (size), protein content, and oil content.

## 1. Significant SNP Summary

- **100SW (size):** 15 genome-wide significant SNPs
- **Protein:** 3 genome-wide significant SNPs
- **Oil:** 8 genome-wide significant SNPs
- **Total unique SNPs:** 26

- **SNPs significant in multiple traits:** 0

## 2. Locus Classification

**Total loci identified:** 24

**Classification breakdown:**
- size_locus: 14 loci (58.3%)
- oil_locus: 7 loci (29.2%)
- protein_locus: 3 loci (12.5%)

**Top loci by number of SNPs:**

| Locus ID | Chromosome | SNPs | Primary Trait | Classification |
|----------|------------|------|---------------|----------------|
| locus_019 | 13 | 2 | oil | oil_locus |
| locus_008 | 15 | 2 | size | size_locus |
| locus_009 | 1 | 1 | size | size_locus |
| locus_003 | 10 | 1 | oil | oil_locus |
| locus_007 | 11 | 1 | protein | protein_locus |

## 3. Cross-Trait Effects

**Effect direction correlations:**
- oil_vs_protein: r = 0.450
- oil_vs_size: r = -0.400
- protein_vs_size: r = -0.643

**Effect direction concordance:**
- oil_vs_protein: 24/26 SNPs (92.3%) have same effect direction
- oil_vs_size: 9/26 SNPs (34.6%) have same effect direction
- protein_vs_size: 7/26 SNPs (26.9%) have same effect direction

## 4. Bayesian Colocalization Results

**Total coloc tests:** 51
**Colocalized loci (PP4 ≥ 0.8):** 0 (0.0%)

**Note:** No loci showed strong evidence of colocalization (PP4 ≥ 0.8).
The maximum PP4 observed was 0.008361.

## 5. Generated Figures

The analysis generated the following visualization files:

1. **cross_trait_scatter.png** - Scatter plots of beta effects between trait pairs
2. **beta_heatmap.png** - Heatmap of standardized beta effects across SNPs
3. **locus_category_barplot.png** - Distribution of locus classifications
4. **coloc_summary_plot.png** - Summary of colocalization results

## 6. Files Generated

| File | Description |
|------|-------------|
| `locus_classification.tsv` | Locus classification results |
| `final_summary_table.tsv` | Combined SNP, locus, and coloc data |
| `figures/` | Directory containing all visualization files |
| `final_analysis_report.md` | This summary report |

## 7. Conclusions

### Key Findings:

1. **Trait-specific genetic architecture:** Most genome-wide significant associations 
   are specific to individual traits with limited pleiotropy.

2. **Limited colocalization:** No strong evidence for shared causal variants 
   between oil, protein, and seed size traits at genome-wide significant loci.

3. **Distinct genetic regulation:** The lack of colocalization suggests 
   independent genetic mechanisms underlying these important agronomic traits.

4. **Methodological framework established:** The complete pipeline from GWAS 
   to colocalization provides a template for future multi-trait analyses.

### Limitations:

1. **Sample size:** Limited power for detecting weak pleiotropic effects
2. **LD considerations:** Colocalization analysis did not incorporate LD structure
3. **Threshold selection:** Genome-wide significance threshold (p < 5e-8) 
   may miss subtler shared genetic effects

### Future Directions:

1. **LD-aware colocalization:** Incorporate LD matrices from genotype data
2. **Fine-mapping:** Identify causal variants within associated loci
3. **Functional validation:** Experimental follow-up of prioritized loci
4. **Multi-omics integration:** Combine with transcriptomics and metabolomics

---

*Report generated on 2026-03-06 21:11:06*
*Analysis pipeline: GWAS Post-Analysis v1.0*