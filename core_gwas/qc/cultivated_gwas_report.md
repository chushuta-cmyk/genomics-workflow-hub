# Cultivated Soybean GWAS Analysis Report

*Generated on: 2026-03-10 17:58:54*

## 1. Overview

This report summarizes the genome-wide association study (GWAS) results for cultivated soybean (Glycine max).

- **Population**: Cultivated soybean
- **Sample size**: 378 accessions
- **Genotyping**: 325,103 SNPs after quality control
- **Traits analyzed**: 100-seed weight (100SW), Protein content, Oil content, log(Oil/Protein) ratio
- **Significance threshold**: p < 5 × 10⁻⁸ (genome-wide)

## 2. Significant SNP Summary

```
Significant SNP Summary
======================
P-value threshold: 5.0e-08
Minimum AF: 0.01
Analysis date: 2026-03-10 16:29:44

Per-trait counts:
----------------------------------------
size           :     15 SNPs
protein        :      3 SNPs
oil            :      8 SNPs
log_ratio      :    126 SNPs

----------------------------------------
Total unique SNPs:    152
Unique SNP IDs:       152
SNPs in >1 trait:       0
```

## 3. Locus Analysis

SNPs within ±250 kb windows were merged into loci.

- **Total loci**: 129
- **SNPs assigned to loci**: 152 (100.0% of significant SNPs)
- **Mean locus size**: 10058 bp
- **Median locus size**: 0 bp
- **Maximum locus size**: 249,232 bp
- **Mean SNPs per locus**: 1.2
- **Multi-trait loci**: 3 (2.3%)

## 4. Linkage Disequilibrium (LD) Analysis

LD decay was calculated for the entire genome (sampling 1 in 1000 SNP pairs).

- **Mean r² at 10,000 bp**: 0.0618
- **Mean r² at 50,000 bp**: 0.0535
- **Mean r² at 100,000 bp**: 0.0549
- **Mean r² at 200,000 bp**: 0.0292
- **Distance where r² < 0.2**: ~0.0 bp

## 5. Haplotype Block Analysis

- **Total haplotype blocks**: 40,902
- **Mean block size**: 5.52 kb
- **Median block size**: 1.21 kb
- **Maximum block size**: 170.60 kb
- **Mean SNPs per block**: 3.17
- **Maximum SNPs per block**: 186

## 6. Files Generated

### GWAS Results
- `gwas/gwas_100SW.tsv` - 100-seed weight GWAS results
- `gwas/gwas_Protein.tsv` - Protein content GWAS results
- `gwas/gwas_Oil.tsv` - Oil content GWAS results
- `gwas/gwas_log_ratio.tsv` - log(Oil/Protein) ratio GWAS results

### Post-GWAS Analysis
- `post_gwas/sig_*.tsv` - Significant SNPs for each trait
- `post_gwas/cross_trait_effect.tsv` - Cross-trait SNP effects
- `post_gwas/locus_summary_*.tsv` - Per-trait locus summaries

### LD Analysis
- `ld_analysis/cultivated_ld_decay_full.ld` - Full genome LD matrix
- `ld_analysis/ld_decay_summary.tsv` - LD decay summary statistics
- `ld_analysis/cultivated_ld_decay_curve.png` - LD decay curve plot
- `ld_analysis/cultivated_haplotype_blocks.blocks` - Haplotype blocks
- `ld_analysis/haplotype_block_summary.tsv` - Block statistics
- `ld_analysis/haplotype_block_*.png` - Block visualization plots

### Visualization
- `reports/cultivated_manhattan_100SW.png` - Manhattan plot for 100SW
- `reports/cultivated_manhattan_Protein.png` - Manhattan plot for Protein
- `reports/cultivated_manhattan_Oil.png` - Manhattan plot for Oil
- `reports/cultivated_manhattan_log_ratio.png` - Manhattan plot for log ratio

## 7. Key Findings

1. **Trait-specific signals**: Each trait showed distinct patterns of genome-wide association.
2. **LD patterns**: Cultivated soybean exhibits characteristic LD decay suitable for GWAS.
3. **Haplotype structure**: The genome is organized into distinct haplotype blocks.
4. **Multi-trait loci**: Several genomic regions contain SNPs associated with multiple traits.

---

*Analysis performed using the cultivated soybean GWAS post-analysis pipeline.*
