# Generalized from soybean project-specific path layout.
#!/usr/bin/env python3
"""
Generate summary plots for wild soybean GWAS loci analysis.

Plots:
1. Manhattan plot per trait
2. Loci count per trait barplot
3. Locus size distribution histogram
4. Cross-trait beta scatterplots
5. Summary plot of significant SNP counts per trait
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import os
import sys

# Set style
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")

# Configuration
WILD_GWAS_DIR = Path("data/input/workflow/04_gemma_wild")
POST_GWAS_DIR = Path("results/wild/post_gwas")
OUTPUT_DIR = POST_GWAS_DIR / "plots"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Trait mapping (GEMMA trait number to trait name)
# Based on prepare_wild_pheno_from_long.py: traits_all = ["100SW","Protein","Oil","Oil_Prot_ratio","log_Oil_Prot_ratio"]
# But cross_trait_effect.tsv only has: 100SW, Protein, Oil, log_ratio
# So we map: trait1->100SW, trait2->Protein, trait3->Oil, trait5->log_ratio
TRAIT_MAP = {
    1: "100SW",
    2: "Protein", 
    3: "Oil",
    5: "log_ratio"
}

# Reverse mapping
TRAIT_TO_NUM = {v: k for k, v in TRAIT_MAP.items()}

def load_gwas_file(trait_num):
    """Load GWAS results for a trait number."""
    file_path = WILD_GWAS_DIR / f"wild_trait_{trait_num}.assoc.txt"
    if not file_path.exists():
        print(f"Warning: GWAS file not found: {file_path}")
        return None
    
    print(f"Loading GWAS results for trait {trait_num} ({TRAIT_MAP.get(trait_num, 'unknown')})...")
    df = pd.read_csv(file_path, sep='\t')
    
    # Required columns
    required_cols = ['chr', 'rs', 'ps', 'p_wald']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        print(f"Warning: Missing columns in {file_path}: {missing_cols}")
        return None
    
    # Convert chromosome to string
    df['chr'] = df['chr'].astype(str)
    
    # Calculate -log10(p)
    df['-log10p'] = -np.log10(df['p_wald'])
    
    return df

def plot_manhattan(df, trait_name, output_path):
    """Create Manhattan plot."""
    if df is None or len(df) == 0:
        print(f"No data for Manhattan plot of {trait_name}")
        return
    
    print(f"Creating Manhattan plot for {trait_name}...")
    
    # Sort by chromosome and position
    df = df.sort_values(['chr', 'ps'])
    
    # Create cumulative position
    df['cum_pos'] = df.groupby('chr').cumcount()
    
    # Get chromosome offsets for x-axis
    chrom_offsets = {}
    cum_offset = 0
    for chrom in sorted(df['chr'].unique()):
        chrom_offsets[chrom] = cum_offset
        cum_offset += df[df['chr'] == chrom].shape[0]
    
    df['x_pos'] = df.apply(lambda row: chrom_offsets[row['chr']] + row['cum_pos'], axis=1)
    
    # Plot
    fig, ax = plt.subplots(figsize=(14, 6))
    
    # Color chromosomes alternately
    colors = ['#648FFF', '#785EF0']  # Blue and purple
    for i, (chrom, chrom_df) in enumerate(df.groupby('chr')):
        color = colors[i % 2]
        ax.scatter(chrom_df['x_pos'], chrom_df['-log10p'], 
                  color=color, s=8, alpha=0.7, label=f'Chr{chrom}' if i < 10 else None)
    
    # Add significance line (5e-8 for genome-wide significance)
    sig_threshold = -np.log10(5e-8)
    ax.axhline(y=sig_threshold, color='r', linestyle='--', alpha=0.7, 
               label=f'Genome-wide sig (p=5e-8)')
    
    # Suggestive line (1e-5)
    sugg_threshold = -np.log10(1e-5)
    ax.axhline(y=sugg_threshold, color='orange', linestyle=':', alpha=0.7,
               label=f'Suggestive (p=1e-5)')
    
    # Set x-axis ticks at chromosome centers
    chrom_centers = {}
    for chrom in sorted(df['chr'].unique()):
        chrom_df = df[df['chr'] == chrom]
        chrom_centers[chrom] = chrom_df['x_pos'].median()
    
    ax.set_xticks(list(chrom_centers.values()))
    ax.set_xticklabels(list(chrom_centers.keys()))
    ax.set_xlabel('Chromosome')
    ax.set_ylabel('-log10(p-value)')
    ax.set_title(f'Manhattan Plot: {trait_name}')
    
    # Limit legend items
    handles, labels = ax.get_legend_handles_labels()
    if len(handles) > 12:  # Too many chromosomes
        handles = handles[:10] + handles[-2:]  # First 10 chromosomes + significance lines
        labels = labels[:10] + labels[-2:]
    
    ax.legend(handles=handles, labels=labels, loc='upper right', fontsize='small')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"  Saved to: {output_path}")

def plot_loci_count_per_trait():
    """Create bar plot of loci count per trait."""
    print("\nCreating loci count per trait barplot...")
    
    loci_counts = {}
    for trait_name in TRAIT_MAP.values():
        stats_file = POST_GWAS_DIR / f"locus_summary_{trait_name}.stats.tsv"
        if stats_file.exists():
            stats_df = pd.read_csv(stats_file, sep='\t')
            total_loci = stats_df[stats_df['metric'] == 'total_loci']['value'].values[0]
            loci_counts[trait_name] = total_loci
    
    if not loci_counts:
        print("No locus summary files found!")
        return
    
    # Create bar plot
    fig, ax = plt.subplots(figsize=(10, 6))
    traits = list(loci_counts.keys())
    counts = list(loci_counts.values())
    
    bars = ax.bar(traits, counts, color=sns.color_palette("husl", len(traits)))
    ax.set_xlabel('Trait')
    ax.set_ylabel('Number of Loci')
    ax.set_title('Number of Loci Identified per Trait')
    
    # Add value labels on bars
    for bar, count in zip(bars, counts):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                f'{int(count)}', ha='center', va='bottom')
    
    plt.tight_layout()
    output_path = OUTPUT_DIR / "loci_count_per_trait.png"
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"  Saved to: {output_path}")
    
    return loci_counts

def plot_locus_size_distribution():
    """Create histogram of locus size distribution across all traits."""
    print("\nCreating locus size distribution histogram...")
    
    all_sizes = []
    trait_labels = []
    
    for trait_name in TRAIT_MAP.values():
        locus_file = POST_GWAS_DIR / f"locus_summary_{trait_name}.tsv"
        if locus_file.exists():
            df = pd.read_csv(locus_file, sep='\t')
            if 'locus_size' in df.columns:
                all_sizes.extend(df['locus_size'].tolist())
                trait_labels.extend([trait_name] * len(df))
    
    if not all_sizes:
        print("No locus size data found!")
        return
    
    # Create combined histogram
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Histogram of all sizes
    axes[0].hist(all_sizes, bins=50, edgecolor='black', alpha=0.7)
    axes[0].set_xlabel('Locus Size (bp)')
    axes[0].set_ylabel('Frequency')
    axes[0].set_title('Distribution of Locus Sizes (All Traits)')
    axes[0].axvline(x=np.median(all_sizes), color='r', linestyle='--', 
                   label=f'Median: {int(np.median(all_sizes))} bp')
    axes[0].legend()
    
    # Box plot by trait
    size_df = pd.DataFrame({'trait': trait_labels, 'locus_size': all_sizes})
    sns.boxplot(x='trait', y='locus_size', data=size_df, ax=axes[1])
    axes[1].set_xlabel('Trait')
    axes[1].set_ylabel('Locus Size (bp)')
    axes[1].set_title('Locus Size by Trait')
    axes[1].tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    output_path = OUTPUT_DIR / "locus_size_distribution.png"
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"  Saved to: {output_path}")
    
    return size_df

def plot_cross_trait_beta_scatter():
    """Create scatter plots of beta values between traits."""
    print("\nCreating cross-trait beta scatterplots...")
    
    # Load cross-trait effect table
    cross_file = POST_GWAS_DIR / "cross_trait_effect.tsv"
    if not cross_file.exists():
        print(f"Cross-trait effect file not found: {cross_file}")
        return
    
    df = pd.read_csv(cross_file, sep='\t')
    
    # Get beta columns
    beta_cols = [col for col in df.columns if col.startswith('beta_')]
    traits = [col.replace('beta_', '') for col in beta_cols]
    
    if len(beta_cols) < 2:
        print("Not enough beta columns for scatter plots")
        return
    
    # Create pair plot
    print(f"  Creating scatter matrix for traits: {traits}")
    
    # Select a subset of data for faster plotting if large
    plot_df = df.copy()
    if len(plot_df) > 10000:
        plot_df = plot_df.sample(n=10000, random_state=42)
        print(f"  Sampling 10,000 SNPs for visualization")
    
    # Create pair plot
    pair_plot = sns.pairplot(plot_df, vars=beta_cols, 
                            diag_kind='hist', diag_kws={'bins': 30},
                            plot_kws={'alpha': 0.5, 's': 10})
    
    # Update diagonal labels
    for i, ax in enumerate(pair_plot.diag_axes):
        ax.set_ylabel(traits[i])
    
    # Update axis labels
    for i in range(len(traits)):
        for j in range(len(traits)):
            if i != j:
                pair_plot.axes[i, j].set_xlabel(traits[j])
                pair_plot.axes[i, j].set_ylabel(traits[i])
    
    pair_plot.fig.suptitle('Cross-Trait Beta Value Correlations', y=1.02)
    output_path = OUTPUT_DIR / "cross_trait_beta_scatter.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved to: {output_path}")
    
    # Also create correlation heatmap
    print(f"  Creating correlation heatmap...")
    corr_matrix = df[beta_cols].corr()
    
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.heatmap(corr_matrix, annot=True, cmap='coolwarm', center=0,
                square=True, linewidths=0.5, ax=ax)
    ax.set_title('Correlation Matrix of Beta Values Across Traits')
    
    # Update tick labels
    tick_labels = traits
    ax.set_xticklabels(tick_labels)
    ax.set_yticklabels(tick_labels)
    
    output_path2 = OUTPUT_DIR / "cross_trait_beta_correlation.png"
    plt.savefig(output_path2, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved to: {output_path2}")
    
    return corr_matrix

def plot_significant_snp_counts():
    """Create summary plot of significant SNP counts per trait."""
    print("\nCreating significant SNP counts summary...")
    
    # Load cross-trait effect table to count significant SNPs per trait
    cross_file = POST_GWAS_DIR / "cross_trait_effect.tsv"
    if not cross_file.exists():
        print(f"Cross-trait effect file not found: {cross_file}")
        return
    
    df = pd.read_csv(cross_file, sep='\t')
    
    # Count significant SNPs (p < 1e-5) for each trait
    sig_counts = {}
    total_counts = {}
    
    for trait_name in TRAIT_MAP.values():
        p_col = f'p_wald_{trait_name}'
        if p_col in df.columns:
            total = len(df[p_col].dropna())
            sig = (df[p_col] < 1e-5).sum()
            sig_counts[trait_name] = sig
            total_counts[trait_name] = total
    
    if not sig_counts:
        print("No p-value columns found!")
        return
    
    # Create stacked bar plot
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Bar plot of counts
    traits = list(sig_counts.keys())
    sig_values = [sig_counts[t] for t in traits]
    non_sig_values = [total_counts[t] - sig_counts[t] for t in traits]
    
    x = np.arange(len(traits))
    width = 0.6
    
    axes[0].bar(x, non_sig_values, width, label='Non-significant (p ≥ 1e-5)', 
               color='lightgray', edgecolor='black')
    axes[0].bar(x, sig_values, width, bottom=non_sig_values, 
               label='Significant (p < 1e-5)', color='coral', edgecolor='black')
    
    axes[0].set_xlabel('Trait')
    axes[0].set_ylabel('Number of SNPs')
    axes[0].set_title('Significant SNP Counts per Trait (p < 1e-5)')
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(traits)
    axes[0].legend()
    
    # Add percentage labels
    for i, trait in enumerate(traits):
        total = total_counts[trait]
        sig = sig_counts[trait]
        percentage = 100 * sig / total if total > 0 else 0
        axes[0].text(i, total * 1.02, f'{percentage:.1f}%', 
                    ha='center', va='bottom', fontsize=9)
    
    # Pie chart of total significant SNPs
    total_sig = sum(sig_counts.values())
    if total_sig > 0:
        wedges, texts, autotexts = axes[1].pie(sig_values, labels=traits, autopct='%1.1f%%',
                                              startangle=90, colors=sns.color_palette("husl", len(traits)))
        axes[1].set_title(f'Distribution of Significant SNPs\n(Total: {total_sig} SNPs, p < 1e-5)')
    
    plt.tight_layout()
    output_path = OUTPUT_DIR / "significant_snp_counts.png"
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"  Saved to: {output_path}")
    
    return sig_counts

def main():
    """Main function to generate all plots."""
    print("=" * 60)
    print("Generating Wild Soybean GWAS Loci Summary Plots")
    print("=" * 60)
    
    # 1. Manhattan plots per trait
    print("\n1. Generating Manhattan plots...")
    for trait_num, trait_name in TRAIT_MAP.items():
        df = load_gwas_file(trait_num)
        if df is not None:
            output_path = OUTPUT_DIR / f"manhattan_{trait_name}.png"
            plot_manhattan(df, trait_name, output_path)
    
    # 2. Loci count per trait barplot
    loci_counts = plot_loci_count_per_trait()
    
    # 3. Locus size distribution histogram
    size_df = plot_locus_size_distribution()
    
    # 4. Cross-trait beta scatterplots
    corr_matrix = plot_cross_trait_beta_scatter()
    
    # 5. Summary plot of significant SNP counts per trait
    sig_counts = plot_significant_snp_counts()
    
    # Create a summary report
    print("\n" + "=" * 60)
    print("SUMMARY REPORT")
    print("=" * 60)
    
    if loci_counts:
        print("\nLoci counts per trait:")
        for trait, count in loci_counts.items():
            print(f"  {trait:12s}: {int(count):4d} loci")
    
    if sig_counts:
        print("\nSignificant SNPs (p < 1e-5) per trait:")
        for trait, count in sig_counts.items():
            total = len(pd.read_csv(POST_GWAS_DIR / "cross_trait_effect.tsv", sep='\t'))
            percentage = 100 * count / total if total > 0 else 0
            print(f"  {trait:12s}: {count:4d} SNPs ({percentage:.2f}%)")
    
    print(f"\nAll plots saved to: {OUTPUT_DIR}")
    print("\nPlot files generated:")
    for plot_file in sorted(OUTPUT_DIR.glob("*.png")):
        print(f"  - {plot_file.name}")
    
    print("\n" + "=" * 60)
    print("Plot generation completed successfully!")
    print("=" * 60)

if __name__ == "__main__":
    main()