#!/usr/bin/env python3
"""
Compare wild vs cultivated GWAS results.
Generate summary tables and plots.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

# Paths
wild_dir = "/data03/karama/projects/soybean_analysis/wild/post_gwas"
cult_dir = "/data03/karama/projects/soybean_analysis/cultivated/post_gwas"
output_dir = "/data03/karama/projects/soybean_analysis/comparison"

# Trait mapping
trait_map = {
    "100SW": "size",
    "Protein": "protein",
    "Oil": "oil",
    "log_ratio": "log_ratio"
}

def load_sig_summary(filepath):
    """Parse sig_summary.txt to get per-trait counts."""
    counts = {}
    with open(filepath, 'r') as f:
        lines = f.readlines()
    in_section = False
    for line in lines:
        line = line.strip()
        if line.startswith("Per-trait counts:"):
            in_section = True
            continue
        if in_section and line.startswith("---"):
            break
        if in_section and ':' in line:
            parts = line.split(':')
            if len(parts) >= 2:
                trait = parts[0].strip()
                # extract number before "SNPs"
                val = parts[1].strip().split()[0]
                try:
                    counts[trait] = int(val)
                except ValueError:
                    pass
    return counts

def load_locus_summary(trait, pop='wild'):
    """Load locus summary file for given trait."""
    if pop == 'wild':
        base = wild_dir
        trait_name = trait
    else:
        base = cult_dir
        # map trait name
        trait_name = trait_map.get(trait, trait.lower())
    # file naming pattern: locus_summary_<trait>.tsv
    # Need to check if file exists with exact name.
    # For cultivated, the file is locus_summary_100SW.tsv even for 'size' trait.
    # Actually the file names are same as wild? Let's try.
    # We'll try both possible names.
    possible_names = [
        f"locus_summary_{trait_name}.tsv",
        f"locus_summary_{trait}.tsv"
    ]
    filepath = None
    for name in possible_names:
        path = os.path.join(base, name)
        if os.path.exists(path):
            filepath = path
            break
    if filepath is None:
        # maybe in significant_snps subdirectory?
        for name in possible_names:
            path = os.path.join(base, "significant_snps", name)
            if os.path.exists(path):
                filepath = path
                break
    if filepath is None:
        print(f"WARNING: Locus summary not found for {trait} {pop}")
        return pd.DataFrame()
    df = pd.read_csv(filepath, sep='\t')
    # ensure chromosome column is string
    df['chr'] = df['chr'].astype(str)
    return df

def loci_overlap(locus1, locus2, max_dist=250000):
    """Check if two loci overlap or are within max_dist."""
    # locus1 and locus2 are rows with chr, start, end, lead_snp_pos
    if locus1['chr'] != locus2['chr']:
        return False
    # if both have start == end (single SNP)
    if locus1['start'] == locus1['end'] and locus2['start'] == locus2['end']:
        # check distance between lead SNPs
        dist = abs(locus1['lead_snp_pos'] - locus2['lead_snp_pos'])
        return dist <= max_dist
    # else check interval overlap
    # extend intervals by max_dist/2? We'll use simple overlap with max_dist padding
    # We'll consider overlap if intervals expanded by max_dist intersect
    s1 = locus1['start'] - max_dist
    e1 = locus1['end'] + max_dist
    s2 = locus2['start'] - max_dist
    e2 = locus2['end'] + max_dist
    return not (e1 < s2 or e2 < s1)

def compare_loci(wild_df, cult_df, trait):
    """Find shared and unique loci between wild and cultivated."""
    if wild_df.empty or cult_df.empty:
        return {'shared': [], 'wild_unique': list(wild_df.index), 'cult_unique': list(cult_df.index)}
    shared = []
    wild_used = set()
    cult_used = set()
    # naive O(N^2) but small N
    for i, wrow in wild_df.iterrows():
        for j, crow in cult_df.iterrows():
            if loci_overlap(wrow, crow):
                shared.append((i, j))
                wild_used.add(i)
                cult_used.add(j)
                break
    wild_unique = [i for i in wild_df.index if i not in wild_used]
    cult_unique = [j for j in cult_df.index if j not in cult_used]
    return {'shared': shared, 'wild_unique': wild_unique, 'cult_unique': cult_unique}

def main():
    # Load significant SNP counts
    wild_sig_counts = load_sig_summary(os.path.join(wild_dir, "significant_snps", "sig_summary.txt"))
    cult_sig_counts = load_sig_summary(os.path.join(cult_dir, "sig_summary.txt"))
    
    # Create summary table
    rows = []
    for trait in trait_map.keys():
        w_count = wild_sig_counts.get(trait, 0)
        c_trait = trait_map[trait]
        c_count = cult_sig_counts.get(c_trait, 0)
        rows.append({
            'trait': trait,
            'wild_sig_snps': w_count,
            'cultivated_sig_snps': c_count,
            'total_unique': w_count + c_count  # not subtracting shared SNPs across populations, just sum
        })
    summary_df = pd.DataFrame(rows)
    summary_path = os.path.join(output_dir, "gwas_comparison_summary.tsv")
    summary_df.to_csv(summary_path, sep='\t', index=False)
    print(f"Saved GWAS comparison summary to {summary_path}")
    
    # Compare loci per trait
    shared_rows = []
    for trait in trait_map.keys():
        wild_loci = load_locus_summary(trait, 'wild')
        cult_loci = load_locus_summary(trait, 'cultivated')
        if wild_loci.empty and cult_loci.empty:
            continue
        result = compare_loci(wild_loci, cult_loci, trait)
        shared = len(result['shared'])
        wild_uniq = len(result['wild_unique'])
        cult_uniq = len(result['cult_unique'])
        shared_rows.append({
            'trait': trait,
            'wild_loci': len(wild_loci),
            'cultivated_loci': len(cult_loci),
            'shared_loci': shared,
            'wild_unique_loci': wild_uniq,
            'cultivated_unique_loci': cult_uniq
        })
        # For each shared locus, record details
        for w_idx, c_idx in result['shared']:
            w = wild_loci.loc[w_idx]
            c = cult_loci.loc[c_idx]
            shared_rows.append({
                'trait': trait,
                'type': 'shared_pair',
                'wild_locus_id': w['locus_id'],
                'wild_chr': w['chr'],
                'wild_start': w['start'],
                'wild_end': w['end'],
                'wild_lead_snp': w['lead_snp'],
                'cultivated_locus_id': c['locus_id'],
                'cultivated_chr': c['chr'],
                'cultivated_start': c['start'],
                'cultivated_end': c['end'],
                'cultivated_lead_snp': c['lead_snp'],
                'distance': abs(w['lead_snp_pos'] - c['lead_snp_pos'])
            })
    shared_df = pd.DataFrame(shared_rows)
    shared_path = os.path.join(output_dir, "shared_loci_summary.tsv")
    shared_df.to_csv(shared_path, sep='\t', index=False)
    print(f"Saved shared loci summary to {shared_path}")
    
    # Plot significant SNP count comparison
    fig, ax = plt.subplots(figsize=(8, 6))
    x = np.arange(len(summary_df))
    width = 0.35
    ax.bar(x - width/2, summary_df['wild_sig_snps'], width, label='Wild', color='skyblue')
    ax.bar(x + width/2, summary_df['cultivated_sig_snps'], width, label='Cultivated', color='lightcoral')
    ax.set_xlabel('Trait')
    ax.set_ylabel('Number of significant SNPs')
    ax.set_title('Significant SNP count comparison (p < 5e-8)')
    ax.set_xticks(x)
    ax.set_xticklabels(summary_df['trait'])
    ax.legend()
    plt.tight_layout()
    count_plot_path = os.path.join(output_dir, "significant_snp_count_comparison.png")
    fig.savefig(count_plot_path, dpi=300)
    print(f"Saved SNP count plot to {count_plot_path}")
    
    # Plot shared loci summary (grouped bar plot)
    fig2, ax2 = plt.subplots(figsize=(10, 6))
    # Filter rows that are not shared_pair
    locus_summary = shared_df[shared_df['type'].isna()] if 'type' in shared_df.columns else shared_df
    if not locus_summary.empty:
        x2 = np.arange(len(locus_summary))
        width2 = 0.2
        ax2.bar(x2 - width2*2, locus_summary['wild_loci'], width2, label='Wild total', color='skyblue')
        ax2.bar(x2 - width2, locus_summary['cultivated_loci'], width2, label='Cultivated total', color='lightcoral')
        ax2.bar(x2, locus_summary['shared_loci'], width2, label='Shared', color='gold')
        ax2.bar(x2 + width2, locus_summary['wild_unique_loci'], width2, label='Wild unique', color='darkblue')
        ax2.bar(x2 + width2*2, locus_summary['cultivated_unique_loci'], width2, label='Cultivated unique', color='darkred')
        ax2.set_xlabel('Trait')
        ax2.set_ylabel('Number of loci')
        ax2.set_title('Locus count comparison')
        ax2.set_xticks(x2)
        ax2.set_xticklabels(locus_summary['trait'])
        ax2.legend()
        plt.tight_layout()
        shared_plot_path = os.path.join(output_dir, "shared_loci_venn_or_barplot.png")
        fig2.savefig(shared_plot_path, dpi=300)
        print(f"Saved shared loci bar plot to {shared_plot_path}")
    else:
        print("No locus summary data to plot.")
    
    # LD decay comparison (Task 1B)
    # Load LD decay summaries
    wild_ld = os.path.join("/data03/karama/projects/soybean_analysis/wild/ld_analysis", "ld_decay_summary.tsv")
    cult_ld = os.path.join("/data03/karama/projects/soybean_analysis/cultivated/ld_analysis", "ld_decay_summary.tsv")
    if os.path.exists(wild_ld) and os.path.exists(cult_ld):
        df_w = pd.read_csv(wild_ld, sep='\t')
        df_c = pd.read_csv(cult_ld, sep='\t')
        # Use distance_bin_start as distance
        df_w['distance'] = df_w['distance_bin_start']
        df_c['distance'] = df_c['distance_bin_start']
        # Plot both curves
        fig3, ax3 = plt.subplots(figsize=(8, 6))
        ax3.plot(df_w['distance'], df_w['mean_r2'], label='Wild', color='blue', marker='o', markersize=3)
        ax3.plot(df_c['distance'], df_c['mean_r2'], label='Cultivated', color='red', marker='s', markersize=3)
        ax3.set_xlabel('Distance (bp)')
        ax3.set_ylabel('Mean r²')
        ax3.set_title('LD decay comparison')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
        # Set x-axis scale maybe log
        ax3.set_xscale('log')
        # Summarize mean r² at specific distances
        distances = [10000, 50000, 100000, 200000, 500000, 1000000]
        summary_rows = []
        for d in distances:
            # find nearest distance in data
            idx_w = (df_w['distance'] - d).abs().argsort()[0]
            idx_c = (df_c['distance'] - d).abs().argsort()[0]
            r2_w = df_w.iloc[idx_w]['mean_r2']
            r2_c = df_c.iloc[idx_c]['mean_r2']
            summary_rows.append({
                'distance_bp': d,
                'wild_mean_r2': r2_w,
                'cultivated_mean_r2': r2_c,
                'difference': r2_w - r2_c
            })
        ld_summary_df = pd.DataFrame(summary_rows)
        ld_summary_path = os.path.join(output_dir, "ld_decay_comparison_summary.tsv")
        ld_summary_df.to_csv(ld_summary_path, sep='\t', index=False)
        print(f"Saved LD decay summary to {ld_summary_path}")
        # Save plot
        ld_plot_path = os.path.join(output_dir, "ld_decay_comparison.png")
        fig3.savefig(ld_plot_path, dpi=300)
        print(f"Saved LD decay comparison plot to {ld_plot_path}")
    else:
        print("LD decay summary files not found. Skipping LD comparison.")
    
    print("Task 1 completed.")

if __name__ == "__main__":
    main()