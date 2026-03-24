# Generalized from soybean project-specific path layout.
#!/usr/bin/env python3
"""
Generate Manhattan plots for GWAS results.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# Set style
try:
    plt.style.use("seaborn")
except Exception:
    pass

def read_gemma_assoc(path):
    """Read GEMMA association results file."""
    df = pd.read_csv(path, delim_whitespace=True, low_memory=False)
    df = df.rename(columns={
        "chr": "CHR",
        "rs": "SNP",
        "ps": "BP",
        "p_wald": "P"
    })
    df = df[df["P"] > 0].copy()
    return df

def manhattan_qq(assoc_file, title="Trait", save_path="manhattan_qq.png"):
    """Generate Manhattan and QQ plots."""
    df = read_gemma_assoc(assoc_file)
    
    # Filter chromosomes 1-20
    df["CHR"] = pd.to_numeric(df["CHR"], errors="coerce")
    df = df[df["CHR"].isin(range(1, 21))].copy()
    df = df.sort_values(["CHR", "BP"])
    
    chr_groups = df.groupby("CHR")
    
    # Calculate cumulative positions
    chr_offsets = {}
    current_offset = 0
    ticks = []
    labels = []
    
    for chrom in range(1, 21):
        if chrom not in chr_groups.groups:
            continue
        group = chr_groups.get_group(chrom)
        chr_offsets[chrom] = current_offset
        max_bp = group["BP"].max()
        ticks.append(current_offset + max_bp / 2)
        labels.append(str(chrom))
        current_offset += max_bp
    
    df["BP_cum"] = df.apply(lambda r: r["BP"] + chr_offsets[r["CHR"]], axis=1)
    df["minus_log10_p"] = -np.log10(df["P"])
    
    # Create figure
    fig = plt.figure(figsize=(10, 8))
    gs = GridSpec(2, 1, height_ratios=[2.5, 1.5])
    
    colors = ["#4C72B0", "#55A868"]
    
    # Manhattan plot
    ax1 = fig.add_subplot(gs[0])
    for i, chrom in enumerate(range(1, 21)):
        if chrom not in chr_groups.groups:
            continue
        group = df[df["CHR"] == chrom]
        ax1.scatter(
            group["BP_cum"],
            group["minus_log10_p"],
            c=colors[i % 2],
            s=6,
            alpha=0.7,
            linewidths=0
        )
    
    # Significance lines
    suggestiveline = 5e-5
    signifline = 5e-8
    ax1.axhline(-np.log10(suggestiveline), color="orange", linestyle="--", linewidth=1)
    ax1.axhline(-np.log10(signifline), color="red", linestyle="-.", linewidth=1)
    
    ax1.set_xticks(ticks)
    ax1.set_xticklabels(labels, fontsize=8)
    ax1.set_ylabel(r"$-log_{10}(P)$")
    ax1.set_title(f"{title} - Manhattan Plot", fontsize=12)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    
    # QQ plot
    ax2 = fig.add_subplot(gs[1])
    p_values = np.sort(df["P"].values)
    n = len(p_values)
    expected = -np.log10((np.arange(1, n + 1) - 0.5) / n)
    observed = -np.log10(p_values)
    
    ax2.scatter(expected, observed, c="#4C72B0", s=8, alpha=0.6, edgecolor="none")
    max_val = max(expected.max(), observed.max())
    ax2.plot([0, max_val], [0, max_val], color="red", linestyle="--", linewidth=1)
    
    ax2.set_xlabel("Expected -log10(P)")
    ax2.set_ylabel("Observed -log10(P)")
    ax2.set_title(f"{title} - QQ Plot", fontsize=12)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    
    fig.subplots_adjust(hspace=0.3)
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {save_path}")

def main():
    """Main function to generate all Manhattan plots."""
    output_dir = "docs/pics"
    os.makedirs(output_dir, exist_ok=True)
    
    # Define traits and corresponding files
    traits = [
        {
            "name": "100SW",
            "file": "data/input/workflow/04_gemma_wild/wild_trait_1.assoc.txt",
            "output": f"{output_dir}/gwas_manhattan_100SW.png"
        },
        {
            "name": "Protein",
            "file": "data/input/workflow/04_gemma_wild/wild_trait_2.assoc.txt",
            "output": f"{output_dir}/gwas_manhattan_Protein.png"
        },
        {
            "name": "Oil",
            "file": "data/input/workflow/04_gemma_wild/wild_trait_3.assoc.txt",
            "output": f"{output_dir}/gwas_manhattan_Oil.png"
        },
        {
            "name": "Oil_Prot_ratio",
            "file": "data/input/workflow/04_gemma_wild/wild_trait_5.assoc.txt",
            "output": f"{output_dir}/gwas_manhattan_Oil_Prot_ratio.png"
        }
    ]
    
    for trait in traits:
        print(f"Generating Manhattan plot for {trait['name']}...")
        if os.path.exists(trait["file"]):
            manhattan_qq(trait["file"], trait["name"], trait["output"])
        else:
            print(f"Warning: File not found: {trait['file']}")

if __name__ == "__main__":
    main()