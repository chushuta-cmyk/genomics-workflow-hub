# Generalized from soybean project-specific path layout.
#!/usr/bin/env python3
import pandas as pd
import numpy as np
import sys

def main():
    blocks_det_file = "results/wild/ld_analysis/wild_haplotype_blocks.blocks.det"
    summary_file = "results/wild/ld_analysis/haplotype_block_summary.tsv"
    
    df = pd.read_csv(blocks_det_file, sep='\s+')
    # Ensure KB column is numeric
    df['KB'] = pd.to_numeric(df['KB'])
    
    total_blocks = len(df)
    mean_size_kb = df['KB'].mean()
    median_size_kb = df['KB'].median()
    max_size_kb = df['KB'].max()
    # Additional stats
    mean_nsnps = df['NSNPS'].mean()
    median_nsnps = df['NSNPS'].median()
    max_nsnps = df['NSNPS'].max()
    
    with open(summary_file, 'w') as fout:
        fout.write("metric\tvalue\n")
        fout.write(f"total_blocks\t{total_blocks}\n")
        fout.write(f"mean_block_size_kb\t{mean_size_kb:.3f}\n")
        fout.write(f"median_block_size_kb\t{median_size_kb:.3f}\n")
        fout.write(f"max_block_size_kb\t{max_size_kb:.3f}\n")
        fout.write(f"mean_snps_per_block\t{mean_nsnps:.2f}\n")
        fout.write(f"median_snps_per_block\t{median_nsnps:.1f}\n")
        fout.write(f"max_snps_per_block\t{max_nsnps}\n")
    
    print(f"Haplotype block summary written to {summary_file}", file=sys.stderr)

if __name__ == "__main__":
    main()