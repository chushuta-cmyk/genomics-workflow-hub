# Generalized from soybean project-specific path layout.
#!/usr/bin/env python3
"""
Generate cultivated soybean GWAS report.
"""

import os
import pandas as pd
import numpy as np
from datetime import datetime

def read_sig_summary(filepath):
    """Read significant SNP summary."""
    with open(filepath, 'r') as f:
        lines = f.readlines()
    return lines

def read_locus_stats(filepath):
    """Read locus statistics."""
    df = pd.read_csv(filepath, sep='\t')
    return df

def read_ld_summary(filepath):
    """Read LD decay summary."""
    df = pd.read_csv(filepath, sep='\t')
    return df

def read_block_summary(filepath):
    """Read haplotype block summary."""
    df = pd.read_csv(filepath, sep='\t')
    return df

def main():
    """Generate comprehensive cultivated soybean GWAS report."""
    # Paths
    base_dir = "results/cultivated"
    sig_summary_file = os.path.join(base_dir, "post_gwas", "sig_summary.txt")
    locus_stats_file = os.path.join(base_dir, "post_gwas", "locus_summary.stats.tsv")
    ld_summary_file = os.path.join(base_dir, "ld_analysis", "ld_decay_summary.tsv")
    block_summary_file = os.path.join(base_dir, "ld_analysis", "haplotype_block_summary.tsv")
    output_file = os.path.join(base_dir, "reports", "cultivated_gwas_report.md")
    
    # Read data
    sig_lines = read_sig_summary(sig_summary_file)
    locus_stats = read_locus_stats(locus_stats_file)
    ld_summary = read_ld_summary(ld_summary_file)
    block_summary = read_block_summary(block_summary_file)
    
    # Convert summary dataframes to dict for easy access
    locus_dict = dict(zip(locus_stats['metric'], locus_stats['value']))
    block_dict = dict(zip(block_summary['metric'], block_summary['value']))
    
    # Generate report
    with open(output_file, 'w') as f:
        f.write("# Cultivated Soybean GWAS Analysis Report\n\n")
        f.write(f"*Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}*\n\n")
        
        f.write("## 1. Overview\n\n")
        f.write("This report summarizes the genome-wide association study (GWAS) results for cultivated soybean (Glycine max).\n\n")
        f.write("- **Population**: Cultivated soybean\n")
        f.write("- **Sample size**: 378 accessions\n")
        f.write("- **Genotyping**: 325,103 SNPs after quality control\n")
        f.write("- **Traits analyzed**: 100-seed weight (100SW), Protein content, Oil content, log(Oil/Protein) ratio\n")
        f.write("- **Significance threshold**: p < 5 × 10⁻⁸ (genome-wide)\n\n")
        
        f.write("## 2. Significant SNP Summary\n\n")
        f.write("```\n")
        for line in sig_lines:
            f.write(line)
        f.write("```\n\n")
        
        f.write("## 3. Locus Analysis\n\n")
        f.write("SNPs within ±250 kb windows were merged into loci.\n\n")
        f.write(f"- **Total loci**: {int(locus_dict.get('total_loci', 0)):,}\n")
        f.write(f"- **SNPs assigned to loci**: {int(locus_dict.get('snps_in_loci', 0)):,} ({locus_dict.get('snps_in_loci', 0)/locus_dict.get('total_snps', 1)*100:.1f}% of significant SNPs)\n")
        f.write(f"- **Mean locus size**: {locus_dict.get('mean_locus_size', 0):.0f} bp\n")
        f.write(f"- **Median locus size**: {locus_dict.get('median_locus_size', 0):.0f} bp\n")
        f.write(f"- **Maximum locus size**: {locus_dict.get('max_locus_size', 0):,.0f} bp\n")
        f.write(f"- **Mean SNPs per locus**: {locus_dict.get('mean_snps_per_locus', 0):.1f}\n")
        f.write(f"- **Multi-trait loci**: {int(locus_dict.get('multi_trait_loci', 0)):,} ({locus_dict.get('multi_trait_loci', 0)/locus_dict.get('total_loci', 1)*100:.1f}%)\n\n")
        
        f.write("## 4. Linkage Disequilibrium (LD) Analysis\n\n")
        f.write("LD decay was calculated for the entire genome (sampling 1 in 1000 SNP pairs).\n\n")
        # Get mean R2 at specific distances
        if len(ld_summary) > 0:
            # Find closest distance bins
            for dist in [10000, 50000, 100000, 200000]:
                # Find bin that contains this distance
                mask = (ld_summary['distance_bin_start'] <= dist) & (ld_summary['distance_bin_end'] >= dist)
                if mask.any():
                    row = ld_summary[mask].iloc[0]
                    if row['mean_r2'] != 'NA':
                        f.write(f"- **Mean r² at {dist:,} bp**: {float(row['mean_r2']):.4f}\n")
            
            # LD decay half-distance (approximate)
            try:
                ld_summary_numeric = ld_summary[ld_summary['mean_r2'] != 'NA'].copy()
                ld_summary_numeric['mean_r2'] = pd.to_numeric(ld_summary_numeric['mean_r2'])
                # Find distance where r2 drops below 0.2
                below_02 = ld_summary_numeric[ld_summary_numeric['mean_r2'] < 0.2]
                if len(below_02) > 0:
                    half_dist = below_02.iloc[0]['distance_bin_start']
                    f.write(f"- **Distance where r² < 0.2**: ~{half_dist:,} bp\n")
            except:
                pass
        
        f.write(f"\n")
        
        f.write("## 5. Haplotype Block Analysis\n\n")
        f.write(f"- **Total haplotype blocks**: {int(block_dict.get('total_blocks', 0)):,}\n")
        f.write(f"- **Mean block size**: {block_dict.get('mean_block_size_kb', 0):.2f} kb\n")
        f.write(f"- **Median block size**: {block_dict.get('median_block_size_kb', 0):.2f} kb\n")
        f.write(f"- **Maximum block size**: {block_dict.get('max_block_size_kb', 0):.2f} kb\n")
        f.write(f"- **Mean SNPs per block**: {block_dict.get('mean_snps_per_block', 0):.2f}\n")
        f.write(f"- **Maximum SNPs per block**: {int(block_dict.get('max_snps_per_block', 0))}\n\n")
        
        f.write("## 6. Files Generated\n\n")
        f.write("### GWAS Results\n")
        f.write("- `gwas/gwas_100SW.tsv` - 100-seed weight GWAS results\n")
        f.write("- `gwas/gwas_Protein.tsv` - Protein content GWAS results\n")
        f.write("- `gwas/gwas_Oil.tsv` - Oil content GWAS results\n")
        f.write("- `gwas/gwas_log_ratio.tsv` - log(Oil/Protein) ratio GWAS results\n\n")
        
        f.write("### Post-GWAS Analysis\n")
        f.write("- `post_gwas/sig_*.tsv` - Significant SNPs for each trait\n")
        f.write("- `post_gwas/cross_trait_effect.tsv` - Cross-trait SNP effects\n")
        f.write("- `post_gwas/locus_summary_*.tsv` - Per-trait locus summaries\n\n")
        
        f.write("### LD Analysis\n")
        f.write("- `ld_analysis/cultivated_ld_decay_full.ld` - Full genome LD matrix\n")
        f.write("- `ld_analysis/ld_decay_summary.tsv` - LD decay summary statistics\n")
        f.write("- `ld_analysis/cultivated_ld_decay_curve.png` - LD decay curve plot\n")
        f.write("- `ld_analysis/cultivated_haplotype_blocks.blocks` - Haplotype blocks\n")
        f.write("- `ld_analysis/haplotype_block_summary.tsv` - Block statistics\n")
        f.write("- `ld_analysis/haplotype_block_*.png` - Block visualization plots\n\n")
        
        f.write("### Visualization\n")
        f.write("- `reports/cultivated_manhattan_100SW.png` - Manhattan plot for 100SW\n")
        f.write("- `reports/cultivated_manhattan_Protein.png` - Manhattan plot for Protein\n")
        f.write("- `reports/cultivated_manhattan_Oil.png` - Manhattan plot for Oil\n")
        f.write("- `reports/cultivated_manhattan_log_ratio.png` - Manhattan plot for log ratio\n\n")
        
        f.write("## 7. Key Findings\n\n")
        f.write("1. **Trait-specific signals**: Each trait showed distinct patterns of genome-wide association.\n")
        f.write("2. **LD patterns**: Cultivated soybean exhibits characteristic LD decay suitable for GWAS.\n")
        f.write("3. **Haplotype structure**: The genome is organized into distinct haplotype blocks.\n")
        f.write("4. **Multi-trait loci**: Several genomic regions contain SNPs associated with multiple traits.\n\n")
        
        f.write("---\n\n")
        f.write("*Analysis performed using the cultivated soybean GWAS post-analysis pipeline.*\n")
    
    print(f"Report generated: {output_file}")

if __name__ == "__main__":
    main()