# Generalized from soybean project-specific path layout.
#!/usr/bin/env python3
"""
Generate MR forest plot for bidirectional analysis.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

def create_forest_plot(mr_results_file, output_file):
    """Create forest plot for bidirectional MR results."""
    
    # Read MR results
    df = pd.read_csv(mr_results_file, sep="\t")
    
    # Filter for key bidirectional comparisons
    # We want: 100SW -> Ratio and Ratio -> 100SW
    # Based on the trace log, we have these relationships
    
    # Extract relevant results
    # From the MR results file, we have:
    # 100SW -> log_ratio (row 6)
    # log_ratio -> 100SW (row 7)
    # 100SW -> Protein (row 4) 
    # 100SW -> Oil (row 2)
    
    # Create a summary DataFrame for the forest plot
    forest_data = []
    
    # Add bidirectional ratio results
    ratio_to_size = df[(df["exposure"] == "log_ratio") & (df["outcome"] == "100SW")]
    size_to_ratio = df[(df["exposure"] == "100SW") & (df["outcome"] == "log_ratio")]
    
    # Add size to protein and oil
    size_to_protein = df[(df["exposure"] == "100SW") & (df["outcome"] == "Protein")]
    size_to_oil = df[(df["exposure"] == "100SW") & (df["outcome"] == "Oil")]
    
    if not ratio_to_size.empty:
        row = ratio_to_size.iloc[0]
        forest_data.append({
            "Comparison": "Ratio → 100SW",
            "Beta": row["beta_ivw"],
            "SE": row["se_ivw"],
            "P": row["p_ivw"],
            "Method": "IVW",
            "Color": "#4C72B0"
        })
    
    if not size_to_ratio.empty:
        row = size_to_ratio.iloc[0]
        forest_data.append({
            "Comparison": "100SW → Ratio",
            "Beta": row["beta_ivw"],
            "SE": row["se_ivw"],
            "P": row["p_ivw"],
            "Method": "IVW",
            "Color": "#55A868"
        })
    
    if not size_to_protein.empty:
        row = size_to_protein.iloc[0]
        forest_data.append({
            "Comparison": "100SW → Protein",
            "Beta": row["beta_ivw"],
            "SE": row["se_ivw"],
            "P": row["p_ivw"],
            "Method": "IVW",
            "Color": "#C44E52"
        })
    
    if not size_to_oil.empty:
        row = size_to_oil.iloc[0]
        forest_data.append({
            "Comparison": "100SW → Oil",
            "Beta": row["beta_ivw"],
            "SE": row["se_ivw"],
            "P": row["p_ivw"],
            "Method": "IVW",
            "Color": "#8172B3"
        })
    
    forest_df = pd.DataFrame(forest_data)
    
    if forest_df.empty:
        print("No MR results found for forest plot.")
        return
    
    # Create forest plot
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Plot each comparison
    for i, row in enumerate(forest_df.itertuples()):
        y_pos = len(forest_df) - i - 1
        beta = row.Beta
        se = row.SE
        
        # Calculate 95% CI
        ci_lower = beta - 1.96 * se
        ci_upper = beta + 1.96 * se
        
        # Plot point and confidence interval
        ax.plot([ci_lower, ci_upper], [y_pos, y_pos], color=row.Color, linewidth=2)
        ax.plot(beta, y_pos, 'o', color=row.Color, markersize=8)
        
        # Add text label with beta and p-value
        p_str = f"{row.P:.2e}" if row.P < 0.001 else f"{row.P:.3f}"
        label = f"β = {beta:.3f} (P = {p_str})"
        ax.text(ci_upper + 0.02 * abs(beta), y_pos, label, 
                va='center', fontsize=9)
    
    # Add vertical line at 0
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
    
    # Set y-axis labels
    ax.set_yticks(range(len(forest_df)))
    ax.set_yticklabels(forest_df["Comparison"].values)
    ax.invert_yaxis()  # Invert so first item is at top
    
    # Set x-axis label
    ax.set_xlabel("Effect Size (β)", fontsize=11)
    
    # Add title
    ax.set_title("Mendelian Randomization - Forest Plot", fontsize=13, pad=15)
    
    # Add grid
    ax.grid(True, alpha=0.3, linestyle='--', axis='x')
    
    # Adjust layout
    plt.tight_layout()
    
    # Save figure
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved forest plot: {output_file}")
    
    # Also print summary
    print("\nMR Results Summary:")
    print("-" * 60)
    for _, row in forest_df.iterrows():
        p_str = f"{row['P']:.2e}" if row['P'] < 0.001 else f"{row['P']:.3f}"
        print(f"{row['Comparison']:20} β = {row['Beta']:8.4f} ± {row['SE']:.4f} (P = {p_str})")

def main():
    """Main function."""
    output_dir = "docs/pics"
    os.makedirs(output_dir, exist_ok=True)
    
    mr_results_file = "data/input/workflow/11_mr_results_wild/mr_ivw_egger_wild_results.tsv"
    output_file = f"{output_dir}/mr_forest_bidirectional.png"
    
    if os.path.exists(mr_results_file):
        create_forest_plot(mr_results_file, output_file)
    else:
        print(f"Error: MR results file not found: {mr_results_file}")

if __name__ == "__main__":
    main()