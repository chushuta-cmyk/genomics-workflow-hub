# Generalized from soybean project-specific path layout.
#!/usr/bin/env python3
"""
Generate MR scatter plots.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def plot_mr_scatter(harmonized_file, mr_results_file, exposure, outcome, output_file):
    """Create MR scatter plot with IVW regression line."""
    
    # Read harmonized data
    if not os.path.exists(harmonized_file):
        print(f"Warning: Harmonized file not found: {harmonized_file}")
        return
    
    df = pd.read_csv(harmonized_file, sep="\t")
    
    # Read MR results to get IVW beta
    mr_df = pd.read_csv(mr_results_file, sep="\t")
    mr_row = mr_df[(mr_df["exposure"] == exposure) & (mr_df["outcome"] == outcome)]
    
    if mr_row.empty:
        print(f"Warning: No MR results for {exposure} -> {outcome}")
        return
    
    beta_ivw = mr_row.iloc[0]["beta_ivw"]
    se_ivw = mr_row.iloc[0]["se_ivw"]
    p_ivw = mr_row.iloc[0]["p_ivw"]
    n_snp = mr_row.iloc[0]["n_snp"]
    
    # Extract data
    bx = df["BETA_exp"].values
    by = df["BETA_out"].values
    
    # Create scatter plot
    plt.figure(figsize=(7, 6))
    
    # Plot scatter points
    plt.scatter(bx, by, alpha=0.5, s=30, edgecolor='none', color='#4C72B0')
    
    # Calculate regression line (IVW through origin)
    x_min, x_max = bx.min(), bx.max()
    x_pad = 0.1 * (x_max - x_min)
    xs = np.linspace(x_min - x_pad, x_max + x_pad, 100)
    ys_ivw = beta_ivw * xs
    
    # Plot IVW regression line
    plt.plot(xs, ys_ivw, color='red', linewidth=2, 
             label=f'IVW: β = {beta_ivw:.4f} (P = {p_ivw:.2e})')
    
    # Add error bars for a few representative points
    # Select 5 representative points for error bars
    if len(df) > 10:
        indices = np.linspace(0, len(df)-1, min(5, len(df)), dtype=int)
        for idx in indices:
            row = df.iloc[idx]
            plt.errorbar(row["BETA_exp"], row["BETA_out"],
                        xerr=row["SE_exp"], yerr=row["SE_out"],
                        color='black', alpha=0.3, linewidth=0.5)
    
    # Format plot
    plt.xlabel(f'SNP effect on {exposure} (β)', fontsize=11)
    plt.ylabel(f'SNP effect on {outcome} (β)', fontsize=11)
    
    # Format title
    p_str = f"{p_ivw:.2e}" if p_ivw < 0.001 else f"{p_ivw:.4f}"
    title = f'MR: {exposure} → {outcome}\nβ = {beta_ivw:.4f} ± {se_ivw:.4f}, P = {p_str}, N SNPs = {n_snp}'
    plt.title(title, fontsize=12, pad=15)
    
    plt.legend(loc='best', fontsize=10)
    plt.grid(True, alpha=0.3, linestyle='--')
    
    # Add zero lines
    plt.axhline(y=0, color='gray', linestyle='-', alpha=0.3, linewidth=0.5)
    plt.axvline(x=0, color='gray', linestyle='-', alpha=0.3, linewidth=0.5)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved scatter plot: {output_file}")

def main():
    """Main function to generate all MR scatter plots."""
    output_dir = "docs/pics"
    os.makedirs(output_dir, exist_ok=True)
    
    mr_results_file = "data/input/workflow/11_mr_results_wild/mr_ivw_egger_wild_results.tsv"
    harm_dir = "data/input/workflow/08_harmonized_wild"
    
    # Define plot pairs based on available data and trace log results
    plot_pairs = [
        # Key bidirectional comparison from trace log
        ("100SW", "log_ratio", "size_to_ratio"),
        ("log_ratio", "100SW", "ratio_to_size"),
        # Other important comparisons
        ("100SW", "Protein", "size_to_protein"),
        ("100SW", "Oil", "size_to_oil"),
    ]
    
    for exposure, outcome, suffix in plot_pairs:
        # Construct harmonized filename
        harm_file = os.path.join(harm_dir, f"{exposure}_to_{outcome}.harmonized.tsv")
        output_file = os.path.join(output_dir, f"mr_scatter_{suffix}.png")
        
        if os.path.exists(harm_file):
            print(f"Generating scatter plot for {exposure} -> {outcome}...")
            plot_mr_scatter(harm_file, mr_results_file, exposure, outcome, output_file)
        else:
            print(f"Warning: Harmonized file not found: {harm_file}")

if __name__ == "__main__":
    main()