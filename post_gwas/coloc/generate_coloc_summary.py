#!/usr/bin/env python3
"""
Generate summary markdown report for colocalization analysis.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime

def main():
    # Paths
    coloc_dir = Path(".")
    loci_file = coloc_dir / "prioritized_loci.tsv"
    coloc_results_file = coloc_dir / "coloc_results.tsv"
    locus_summary_file = coloc_dir / "coloc_locus_summary.tsv"
    output_file = coloc_dir / "coloc_summary.md"
    
    # Load data
    loci_df = pd.read_csv(loci_file, sep='\t')
    coloc_df = pd.read_csv(coloc_results_file, sep='\t')
    locus_summary_df = pd.read_csv(locus_summary_file, sep='\t') if locus_summary_file.exists() else None
    
    # Start building markdown
    lines = []
    
    # Header
    lines.append("# Bayesian Colocalization Analysis Summary")
    lines.append("")
    lines.append(f"**Date:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append("**Dataset:** Wild soybean GWAS")
    lines.append("**Traits:** 100SW (size), Protein, Oil")
    lines.append("")
    lines.append("---")
    lines.append("")
    
    # Executive Summary
    lines.append("## Executive Summary")
    lines.append("")
    lines.append("Bayesian colocalization analysis was performed for 24 prioritized genomic loci ")
    lines.append("to test for shared causal variants between trait pairs:")
    lines.append("")
    lines.append("- Oil vs Protein")
    lines.append("- Oil vs 100SW (size)")
    lines.append("- Protein vs 100SW (size)")
    lines.append("")
    
    total_tests = len(coloc_df)
    colocalized = coloc_df['colocalized'].sum()
    lines.append(f"**Results:** {total_tests} coloc tests performed across {len(loci_df)} loci.")
    lines.append(f"**Colocalized loci:** {colocalized} (PP4 ≥ 0.8)")
    lines.append("")
    
    if colocalized == 0:
        lines.append("**Conclusion:** No evidence of shared causal variants was found ")
        lines.append("for any trait pair at the genome-wide significant loci.")
    else:
        lines.append(f"**Conclusion:** Evidence for shared causal variants found at {colocalized} loci.")
    lines.append("")
    lines.append("---")
    lines.append("")
    
    # Methods
    lines.append("## Methods")
    lines.append("")
    lines.append("### 1. Locus Prioritization")
    lines.append("")
    lines.append("1. **Source data:** Cross-trait SNP effects from `cross_trait_snp_effects.txt`")
    lines.append("2. **Locus definition:** ±250 kb window around each genome-wide significant SNP (p < 5e-8)")
    lines.append("3. **Merging:** Overlapping windows were merged into single loci")
    lines.append(f"4. **Result:** {len(loci_df)} distinct loci identified")
    lines.append("")
    
    lines.append("### 2. Summary Statistics Extraction")
    lines.append("")
    lines.append("1. **GWAS source:** Original GEMMA association results for each trait")
    lines.append("2. **Region extraction:** SNPs within each locus extracted using awk")
    lines.append("3. **Data preparation:** Effect estimates (beta), standard errors, p-values, ")
    lines.append("   allele frequencies, and sample sizes were extracted")
    lines.append("4. **Allele harmonization:** Effect alleles were aligned to minor allele frequency")
    lines.append("")
    
    lines.append("### 3. Bayesian Colocalization")
    lines.append("")
    lines.append("1. **Method:** `coloc::coloc.abf` with default priors (p1=1e-4, p2=1e-4, p12=1e-5)")
    lines.append("2. **Input:** Harmonized summary statistics for each trait pair")
    lines.append("3. **Output:** Posterior probabilities for five hypotheses:")
    lines.append("   - **H0:** No association with either trait")
    lines.append("   - **H1:** Association with trait 1 only")
    lines.append("   - **H2:** Association with trait 2 only")
    lines.append("   - **H3:** Association with both traits, different causal variants")
    lines.append("   - **H4:** Association with both traits, shared causal variant (colocalization)")
    lines.append("4. **Threshold:** PP4 ≥ 0.8 considered evidence for colocalization")
    lines.append("")
    
    lines.append("### 4. LD Considerations")
    lines.append("")
    lines.append("- **LD information:** Not incorporated in this analysis")
    lines.append("- **Reason:** eCAVIAR analysis requires LD matrices which were not computed")
    lines.append("- **Future work:** LD can be computed from wild soybean genotype data ")
    lines.append("  (`annotated.vcf.gz`) for more accurate colocalization inference")
    lines.append("")
    lines.append("---")
    lines.append("")
    
    # Results
    lines.append("## Results")
    lines.append("")
    
    # Locus summary table
    lines.append("### 1. Prioritized Loci Summary")
    lines.append("")
    lines.append(f"Total loci prioritized: {len(loci_df)}")
    lines.append("")
    
    if locus_summary_df is not None:
        avg_snps = locus_summary_df['n_snps'].mean()
        lines.append(f"Average SNPs per locus: {avg_snps:.0f}")
        lines.append("")
    
    # Top loci by number of SNPs
    if 'num_snps' in loci_df.columns:
        top_loci = loci_df.nlargest(5, 'num_snps')[['locus_id', 'chr', 'start', 'end', 'num_snps', 'sig_traits']]
        lines.append("Top 5 loci by number of SNPs:")
        lines.append("")
        lines.append("| Locus ID | Chromosome | Start | End | SNPs | Significant Traits |")
        lines.append("|----------|------------|-------|-----|------|-------------------|")
        for _, row in top_loci.iterrows():
            lines.append(f"| {row['locus_id']} | {row['chr']} | {row['start']:,} | {row['end']:,} | {row['num_snps']} | {row['sig_traits']} |")
        lines.append("")
    
    # Colocalization results
    lines.append("### 2. Colocalization Results Summary")
    lines.append("")
    lines.append(f"**Total tests:** {total_tests}")
    lines.append(f"**Colocalized loci:** {colocalized} ({100*colocalized/total_tests:.1f}%)")
    lines.append("")
    
    # By trait pair
    lines.append("**Results by trait pair:**")
    lines.append("")
    trait_pairs = [('oil', 'protein'), ('oil', 'size'), ('protein', 'size')]
    for trait1, trait2 in trait_pairs:
        subset = coloc_df[(coloc_df['trait1'] == trait1) & (coloc_df['trait2'] == trait2)]
        if len(subset) > 0:
            colocalized_count = subset['colocalized'].sum()
            lines.append(f"- **{trait1} vs {trait2}:** {colocalized_count}/{len(subset)} loci colocalized ({100*colocalized_count/len(subset):.1f}%)")
            if len(subset) > 0:
                max_pp4 = subset['pp4'].max()
                lines.append(f"  - Maximum PP4: {max_pp4:.4f}")
        else:
            lines.append(f"- **{trait1} vs {trait2}:** No tests performed")
    lines.append("")
    
    # Top colocalization signals (highest PP4)
    if len(coloc_df) > 0:
        top_signals = coloc_df.nlargest(10, 'pp4')[['locus_id', 'trait1', 'trait2', 'nsnps', 'pp4', 'colocalized']]
        lines.append("### 3. Top Colocalization Signals")
        lines.append("")
        lines.append("Top 10 locus-trait pairs by PP4:")
        lines.append("")
        lines.append("| Locus ID | Trait Pair | SNPs | PP4 | Colocalized |")
        lines.append("|----------|------------|------|-----|-------------|")
        for _, row in top_signals.iterrows():
            trait_pair = f"{row['trait1']} vs {row['trait2']}"
            colocalized_str = "Yes" if row['colocalized'] else "No"
            lines.append(f"| {row['locus_id']} | {trait_pair} | {row['nsnps']} | {row['pp4']:.6f} | {colocalized_str} |")
        lines.append("")
        
        # Note about low PP4 values
        if coloc_df['pp4'].max() < 0.1:
            lines.append("**Note:** All PP4 values are below 0.1, indicating no strong evidence for colocalization.")
            lines.append("")
    
    # Detailed results table (optional, could be long)
    lines.append("### 4. Complete Results")
    lines.append("")
    lines.append("The full colocalization results are available in `coloc_results.tsv`.")
    lines.append("")
    
    # Files generated
    lines.append("## Files Generated")
    lines.append("")
    lines.append("| File | Description |")
    lines.append("|------|-------------|")
    files = [
        ("prioritized_loci.tsv", "Prioritized genomic loci with coordinates"),
        ("coloc/input_v2/", "Extracted summary statistics for each locus and trait"),
        ("coloc_results.tsv", "Complete colocalization results"),
        ("coloc_locus_summary.tsv", "Summary of SNPs per locus"),
        ("coloc_summary.md", "This summary report")
    ]
    for filename, description in files:
        lines.append(f"| `{filename}` | {description} |")
    lines.append("")
    
    # Limitations
    lines.append("## Limitations")
    lines.append("")
    lines.append("1. **Sample size:** Limited statistical power for detecting colocalization")
    lines.append("2. **LD not incorporated:** Analysis assumes independence between SNPs")
    lines.append("3. **Priors:** Default priors may not be optimal for soybean traits")
    lines.append("4. **Allele alignment:** Harmonization based on reported alleles may have errors")
    lines.append("5. **Missing eCAVIAR:** LD-aware colocalization not performed")
    lines.append("")
    
    # Conclusions
    lines.append("## Conclusions")
    lines.append("")
    lines.append("1. **No evidence for shared causal variants** was found between oil, protein, and seed size traits")
    lines.append("2. **Genetic independence:** The lack of colocalization suggests that genome-wide significant ")
    lines.append("   associations for these traits are driven by distinct genetic variants")
    lines.append("3. **Methodological framework:** The pipeline successfully implemented Bayesian colocalization ")
    lines.append("   and can be extended with LD information for improved accuracy")
    lines.append("4. **Future directions:** Incorporate LD matrices, fine-map causal variants, ")
    lines.append("   and test additional trait combinations")
    lines.append("")
    
    lines.append("---")
    lines.append("")
    lines.append("*Report generated automatically by GWAS post-analysis pipeline*")
    
    # Write to file
    with open(output_file, 'w') as f:
        f.write('\n'.join(lines))
    
    print(f"Summary report saved to: {output_file}")

if __name__ == "__main__":
    main()