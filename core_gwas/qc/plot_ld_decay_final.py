# Generalized from soybean project-specific path layout.
#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

def main():
    summary_file = "results/wild/ld_analysis/ld_decay_summary.tsv"
    plot_file = "results/wild/ld_analysis/wild_ld_decay_curve.png"
    
    df = pd.read_csv(summary_file, sep='\t')
    # Filter out NA rows
    df = df[df['mean_r2'] != 'NA']
    df['mean_r2'] = pd.to_numeric(df['mean_r2'])
    # Use mid-point of bin as distance
    df['distance'] = (df['distance_bin_start'] + df['distance_bin_end']) / 2
    
    plt.figure(figsize=(10, 6))
    plt.plot(df['distance'], df['mean_r2'], 'b-', linewidth=2, label='Mean $r^2$')
    plt.xlabel('Distance (bp)', fontsize=12)
    plt.ylabel('Mean $r^2$', fontsize=12)
    plt.title('LD Decay Curve (Wild Soybean)', fontsize=14)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.xlim(0, df['distance'].max())
    # Set y-axis limit with some padding
    ymax = df['mean_r2'].max() * 1.1
    plt.ylim(0, ymax)
    # Add legend
    plt.legend()
    plt.tight_layout()
    plt.savefig(plot_file, dpi=300)
    print(f"LD decay plot saved to {plot_file}", file=sys.stderr)

if __name__ == "__main__":
    main()