# Generalized from soybean project-specific path layout.
#!/usr/bin/env python3
"""
Generate MR scatter plots for cultivated soybean MR results.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import stats

def create_scatter_plot(ratio_size_file, output_dir):
    """Create scatter plots for ratio-size bidirectional MR."""
    
    if not os.path.exists(ratio_size_file):
        print(f"Error: MR effects file not found: {ratio_size_file}")
        return
    
    # Read SNP-level data
    df = pd.read_csv(ratio_size_file, sep="\t")
    
    if df.empty:
        print("No SNP data found.")
        return
    
    # Create scatter plot for ratio vs size
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Plot 1: Ratio vs Size (for 100SW → Ratio direction)
    x1 = df["beta_size"].values
    y1 = df["beta_ratio"].values
    
    # Calculate regression line
    slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(x1, y1)
    x_line1 = np.array([x1.min(), x1.max()])
    y_line1 = slope1 * x_line1 + intercept1
    
    ax1.scatter(x1, y1, alpha=0.6, s=20, color='#4C72B0', edgecolor='black', linewidth=0.5)
    ax1.plot(x_line1, y_line1, 'r--', linewidth=2, label=f'Slope = {slope1:.3f}\nR² = {r_value1**2:.3f}')
    
    ax1.set_xlabel("SNP effect on 100SW", fontsize=12, fontweight='bold')
    ax1.set_ylabel("SNP effect on Oil:Protein Ratio", fontsize=12, fontweight='bold')
    ax1.set_title("Cultivated Soybean: 100SW → Ratio MR Scatter", fontsize=13, fontweight='bold')
    ax1.legend(loc='best')
    ax1.grid(True, alpha=0.3, linestyle='--')
    
    # Plot 2: Size vs Ratio (for Ratio → 100SW direction)
    x2 = df["beta_ratio"].values
    y2 = df["beta_size"].values
    
    # Calculate regression line
    slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(x2, y2)
    x_line2 = np.array([x2.min(), x2.max()])
    y_line2 = slope2 * x_line2 + intercept2
    
    ax2.scatter(x2, y2, alpha=0.6, s=20, color='#55A868', edgecolor='black', linewidth=0.5)
    ax2.plot(x_line2, y_line2, 'r--', linewidth=2, label=f'Slope = {slope2:.3f}\nR² = {r_value2**2:.3f}')
    
    ax2.set_xlabel("SNP effect on Oil:Protein Ratio", fontsize=12, fontweight='bold')
    ax2.set_ylabel("SNP effect on 100SW", fontsize=12, fontweight='bold')
    ax2.set_title("Cultivated Soybean: Ratio → 100SW MR Scatter", fontsize=13, fontweight='bold')
    ax2.legend(loc='best')
    ax2.grid(True, alpha=0.3, linestyle='--')
    
    plt.tight_layout()
    
    # Save combined plot
    output_file = f"{output_dir}/cultivated_mr_scatter_ratio_size.png"
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved scatter plot: {output_file}")
    
    # Create individual scatter plots for each direction
    for direction, x_col, y_col, color, title, suffix in [
        ("100SW_to_Ratio", "beta_size", "beta_ratio", "#4C72B0", 
         "Cultivated Soybean: 100SW → Oil:Protein Ratio MR Scatter", "size_to_ratio"),
        ("Ratio_to_100SW", "beta_ratio", "beta_size", "#55A868",
         "Cultivated Soybean: Oil:Protein Ratio → 100SW MR Scatter", "ratio_to_size")
    ]:
        fig, ax = plt.subplots(figsize=(8, 6))
        x = df[x_col].values
        y = df[y_col].values
        
        # Calculate regression line
        slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
        x_line = np.array([x.min(), x.max()])
        y_line = slope * x_line + intercept
        
        ax.scatter(x, y, alpha=0.6, s=25, color=color, edgecolor='black', linewidth=0.5)
        ax.plot(x_line, y_line, 'r--', linewidth=2, 
                label=f'Slope = {slope:.3f}\nR² = {r_value**2:.3f}\nP = {p_value:.2e}')
        
        ax.set_xlabel(f"SNP effect on {x_col.replace('beta_', '')}", fontsize=12, fontweight='bold')
        ax.set_ylabel(f"SNP effect on {y_col.replace('beta_', '')}", fontsize=12, fontweight='bold')
        ax.set_title(title, fontsize=13, fontweight='bold')
        ax.legend(loc='best')
        ax.grid(True, alpha=0.3, linestyle='--')
        
        plt.tight_layout()
        output_file = f"{output_dir}/cultivated_mr_scatter_{suffix}.png"
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved scatter plot: {output_file}")

def main():
    """Main function."""
    output_dir = "docs/pics"
    os.makedirs(output_dir, exist_ok=True)
    
    # MR effects file for ratio-size
    ratio_size_file = "data/input/workflow/output/04_gemma_results/MR_ratio_size_effects.txt"
    
    create_scatter_plot(ratio_size_file, output_dir)

if __name__ == "__main__":
    main()