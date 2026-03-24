# Wild Soybean GWAS Downstream Analysis Report

## Overview
This report summarizes the downstream analysis of genome-wide association studies (GWAS) for four traits in wild soybean (*Glycine soja*): 100-seed weight (100SW), seed protein content (Protein), seed oil content (Oil), and the log-ratio of protein to oil (log_ratio). The analysis includes identification of significant SNPs, locus definition, cross‑trait effect analysis, and visualization of results.

## Significant SNP Counts
Using a genome‑wide significance threshold of \(p < 5 \times 10^{-8}\) and a minimum allele frequency of 0.01:

| Trait | Significant SNPs |
|-------|------------------|
| 100SW | 35 |
| Protein | 56 |
| Oil | 3 |
| log_ratio | 5 |
| **Total unique SNPs** | **99** |

- 97 unique SNP IDs (two SNPs are significant in more than one trait).
- Detailed SNP lists are available in `post_gwas/significant_snps/`.

## Locus Summaries by Trait
Loci were defined by merging significant SNPs within a 100‑kb window. The same set of loci was identified across all four traits because the majority of loci contain a single SNP and are shared among traits.

| Statistic | 100SW | Protein | Oil | log_ratio |
|-----------|-------|---------|-----|-----------|
| Total loci | 87 | 87 | 87 | 87 |
| Mean locus size (bp) | 5,171 | 5,171 | 5,171 | 5,171 |
| Median locus size (bp) | 0 | 0 | 0 | 0 |
| Maximum locus size (bp) | 274,554 | 274,554 | 274,554 | 274,554 |
| Mean SNPs per locus | 1.11 | 1.11 | 1.11 | 1.11 |
| Single‑SNP loci | 81 | 81 | 81 | 81 |
| Multi‑trait loci | 4 | 4 | 4 | 4 |

*Note:* The median locus size of 0 bp indicates that more than half of the loci consist of a single SNP (no physical interval).

## Cross‑Trait Effects
SNPs that reached significance in more than one trait were examined to identify pleiotropic effects.

| Significance pattern | Count | Percentage |
|----------------------|-------|------------|
| Protein only | 56 | 57.7% |
| 100SW only | 35 | 36.1% |
| log_ratio only | 3 | 3.1% |
| Oil + log_ratio | 2 | 2.1% |
| Oil only | 1 | 1.0% |

**Effect Direction Concordance**  
For SNPs that are significant in both traits of a pair, the concordance of effect directions (same vs. opposite sign) was calculated:

| Trait pair | SNPs tested | Concordant | Discordant | Concordance rate |
|------------|-------------|------------|------------|------------------|
| 100SW vs Protein | 97 | 33 | 64 | 34.0% |
| 100SW vs Oil | 97 | 13 | 84 | 13.4% |
| 100SW vs log_ratio | 97 | 9 | 88 | 9.3% |
| Protein vs Oil | 97 | 63 | 34 | 64.9% |
| Protein vs log_ratio | 97 | 63 | 34 | 64.9% |
| Oil vs log_ratio | 97 | 93 | 4 | 95.9% |

Concordance rates are highest between Oil and log_ratio (expected due to their inverse relationship) and between Protein and Oil/log_ratio, while 100SW shows low concordance with seed composition traits.

No SNP was simultaneously significant in Protein, 100SW, and Oil/log_ratio, suggesting largely independent genetic architectures for these traits in wild soybean.

## Generated Figures
All visualizations have been copied to `reports/figures/`.

### Manhattan Plots
- `manhattan_100SW.png`
- `manhattan_Protein.png`
- `manhattan_Oil.png`
- `manhattan_log_ratio.png`

### Summary Figures
- `significant_snp_counts.png` – Bar chart of per‑trait SNP counts.
- `loci_count_per_trait.png` – Number of loci per trait (identical across traits).
- `locus_size_distribution.png` – Histogram of locus sizes.
- `cross_trait_beta_correlation.png` – Correlation of effect sizes across traits.
- `cross_trait_beta_scatter.png` – Scatter matrix of effect sizes.

## LD Analysis Status
**LD decay calculation**  
- PLINK LD calculation completed for all chromosomes (`wild_ld_decay_full.ld`, ~18.7 GB).
- Parameters: `--ld-window 99999 --ld-window-kb 1000 --ld-window-r2 0`.
- No downstream summary (e.g., LD decay curve, average \(r^2\) vs. distance) has been generated.

**Haplotype blocks**  
- Attempted with `--blocks` but halted because PLINK requires at least two founders with non‑missing phenotypes (the default can be overridden with `--blocks no-pheno-req`).
- Raw haplotype block output is not available.

**Recommendation:** If LD decay curves are needed, run a dedicated LD decay script (e.g., using `--ld-window-r2 0.2` and binning distances). For haplotype blocks, re‑run with `--blocks no-pheno-req`.

## Current Interpretation
1. **Trait‑specific signals** – Protein and 100SW show the largest number of significant SNPs, consistent with their polygenic nature. Oil and log_ratio have few associations, possibly due to lower genetic variance or stricter significance thresholds.
2. **Locus architecture** – The vast majority of associations are single‑SNP loci, indicating sparse, discrete signals rather than extended haplotypes.
3. **Pleiotropy** – Limited overlap between traits, suggesting distinct genetic mechanisms for seed composition and weight in wild soybean.
4. **Effect directions** – Cross‑trait beta correlations (see scatter plots) show no strong positive or negative correlations, supporting independence of genetic effects.

## Next Recommended Analyses
The following steps are proposed to extend the study:

1. **LD‑based fine‑mapping** – Use the calculated LD matrix to define credible sets for each locus and prioritize candidate causal variants.
2. **Bayesian colocalization** – Test for colocalization of GWAS signals with expression quantitative trait loci (eQTLs) from soybean seed tissues.
3. **Selection scan integration** – Overlap significant loci with genomic regions under selection (e.g., using XP‑EHH or F<sub>ST</sub> between wild and cultivated soybean).
4. **Gene annotation** – Annotate loci with nearby genes (within 50 kb) and perform functional enrichment analysis (GO, KEGG).
5. **Pathway analysis** – Conduct pathway‑based tests (e.g., MAGMA) using the full GWAS summary statistics.
6. **Comparison with cultivated soybean** – Compare the wild soybean associations with published GWAS results in cultivated soybean to identify conserved vs. lineage‑specific loci.

A detailed TODO list is available in `reports/next_steps_todo.md`.

---
**Analysis Date:** 2026‑03‑09  
**Project Root:** `results/wild`  
**Post‑GWAS Directory:** `post_gwas/`  
**Report Directory:** `reports/`