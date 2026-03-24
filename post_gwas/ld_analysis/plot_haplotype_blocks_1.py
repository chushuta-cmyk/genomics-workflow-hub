# Generalized from soybean project-specific path layout.
#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

def main():
    blocks_det_file = "results/cultivated/ld_analysis/cultivated_haplotype_blocks.blocks.det"
    
    df = pd.read_csv(blocks_det_file, sep='\s+')
    df['KB'] = pd.to_numeric(df['KB'])
    
    # Plot 1: Block size distribution (histogram)
    plt.figure(figsize=(10, 6))
    plt.hist(df['KB'], bins=50, edgecolor='black', alpha=0.7)
    plt.xlabel('Block Size (KB)', fontsize=12)
    plt.ylabel('Frequency', fontsize=12)
    plt.title('Haplotype Block Size Distribution (Cultivated Soybean)', fontsize=14)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig('results/cultivated/ld_analysis/haplotype_block_size_distribution.png', dpi=300)
    print("Saved block size distribution plot", file=sys.stderr)
    
    # Plot 2: Block count per chromosome
    plt.figure(figsize=(12, 6))
    block_counts = df['CHR'].value_counts().sort_index()
    # Ensure chromosome order as integer if possible
    try:
        block_counts.index = block_counts.index.astype(int)
        block_counts = block_counts.sort_index()
    except:
        pass
    bars = plt.bar(block_counts.index.astype(str), block_counts.values, color='steelblue')
    plt.xlabel('Chromosome', fontsize=12)
    plt.ylabel('Number of Haplotype Blocks', fontsize=12)
    plt.title('Haplotype Block Count per Chromosome (Cultivated Soybean)', fontsize=14)
    # Add count labels on bars
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                 f'{int(height)}', ha='center', va='bottom', fontsize=9)
    plt.grid(True, linestyle='--', alpha=0.5, axis='y')
    plt.tight_layout()
    plt.savefig('results/cultivated/ld_analysis/haplotype_block_count_per_chr.png', dpi=300)
    print("Saved block count per chromosome plot", file=sys.stderr)

if __name__ == "__main__":
    main()