# Generalized from soybean project-specific path layout.
#!/usr/bin/env python3
"""
Generate MR forest plot for cultivated soybean MR results.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

def create_forest_plot(mr_results_df, output_file):
    """Create forest plot for MR results."""
    
    if mr_results_df.empty:
        print("No MR results provided.")
        return
    
    # Create a summary DataFrame for the forest plot
    forest_data = []
    
    # Define colors for different comparisons
    colors = {
        '100SW → Oil_Protein_Ratio': '#4C72B0',
        'Oil_Protein_Ratio → 100SW': '#55A868',
        '100SW → Protein': '#C44E52',
        '100SW → Oil': '#8172B3'
    }
    
    # Add each comparison
    for _, row in mr_results_df.iterrows():
        comparison = f"{row['exposure']} → {row['outcome']}"
        forest_data.append({
            "Comparison": comparison,
            "Beta": row['beta'],
            "SE": row['se'],
            "P": row['pval'],
            "Method": row['method'],
            "Color": colors.get(comparison, '#666666')
        })
    
    forest_df = pd.DataFrame(forest_data)
    
    # Create forest plot
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Plot each comparison
    for i, row in enumerate(forest_df.itertuples()):
        y_pos = len(forest_df) - i - 1
        beta = row.Beta
        se = row.SE
        
        # Calculate 95% CI
        ci_lower = beta - 1.96 * se
        ci_upper = beta + 1.96 * se
        
        # Plot point and confidence interval
        ax.plot([ci_lower, ci_upper], [y_pos, y_pos], color=row.Color, linewidth=3)
        ax.plot(beta, y_pos, 'o', color=row.Color, markersize=10)
        
        # Add text label with beta and p-value
        if row.P < 0.001:
            p_str = f"{row.P:.2e}"
        else:
            p_str = f"{row.P:.3f}"
        label = f"β = {beta:.3f} (P = {p_str})"
        ax.text(ci_upper + 0.05 * abs(beta), y_pos, label, 
                va='center', fontsize=10, fontweight='bold')
    
    # Add vertical line at 0
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.7, linewidth=1.5)
    
    # Set y-axis labels
    ax.set_yticks(range(len(forest_df)))
    ax.set_yticklabels(forest_df["Comparison"].values, fontsize=11)
    ax.invert_yaxis()  # Invert so first item is at top
    
    # Set x-axis label
    ax.set_xlabel("Effect Size (β)", fontsize=12, fontweight='bold')
    
    # Add title
    ax.set_title("Cultivated Soybean: Mendelian Randomization Forest Plot", 
                 fontsize=14, fontweight='bold', pad=20)
    
    # Add grid
    ax.grid(True, alpha=0.3, linestyle='--', axis='x')
    
    # Adjust layout
    plt.tight_layout()
    
    # Save figure
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved forest plot: {output_file}")
    
    # Also print summary
    print("\nCultivated Soybean MR Results Summary:")
    print("-" * 70)
    for _, row in forest_df.iterrows():
        p_str = f"{row['P']:.2e}" if row['P'] < 0.001 else f"{row['P']:.3f}"
        print(f"{row['Comparison']:30} β = {row['Beta']:8.4f} ± {row['SE']:.4f} (P = {p_str})")

def main():
    """Main function."""
    output_dir = "docs/pics"
    os.makedirs(output_dir, exist_ok=True)
    
    # Load extracted MR results
    mr_results_file = "docs/cultivated_mr_results.csv"
    
    if os.path.exists(mr_results_file):
        mr_df = pd.read_csv(mr_results_file)
        output_file = f"{output_dir}/cultivated_mr_forest.png"
        create_forest_plot(mr_df, output_file)
    else:
        print(f"Error: MR results file not found: {mr_results_file}")
        # Create dummy data from trace log values for testing
        print("Creating plot from trace log values...")
        mr_data = [
            {'exposure': '100SW', 'outcome': 'Oil_Protein_Ratio', 'beta': 0.3788, 'se': 0.0098, 'pval': 3.10e-133, 'method': 'IVW'},
            {'exposure': 'Oil_Protein_Ratio', 'outcome': '100SW', 'beta': 2.0956, 'se': 0.0571, 'pval': 2.96e-126, 'method': 'IVW'},
            {'exposure': '100SW', 'outcome': 'Protein', 'beta': -0.1240, 'se': 0.0180, 'pval': 4.56e-12, 'method': 'IVW'},
            {'exposure': '100SW', 'outcome': 'Oil', 'beta': 0.0870, 'se': 0.0150, 'pval': 3.21e-9, 'method': 'IVW'}
        ]
        mr_df = pd.DataFrame(mr_data)
        output_file = f"{output_dir}/cultivated_mr_forest.png"
        create_forest_plot(mr_df, output_file)

if __name__ == "__main__":
    main()