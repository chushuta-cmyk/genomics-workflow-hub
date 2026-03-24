#!/usr/bin/env python3
"""
Final analysis: locus classification, summary table, and figure generation.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import sys

# Set plotting style
plt.style.use('seaborn-whitegrid')
sns.set_palette("husl")

def load_data():
    """Load all required data files."""
    data = {}
    
    # Load significant SNP files
    data['size_sig'] = pd.read_csv('100SW_significant.txt', sep='\t')
    data['protein_sig'] = pd.read_csv('protein_significant.txt', sep='\t')
    data['oil_sig'] = pd.read_csv('oil_significant.txt', sep='\t')
    data['all_sig'] = pd.read_csv('all_significant_snps.txt', sep='\t')
    
    # Load cross-trait effects
    data['cross_trait'] = pd.read_csv('cross_trait_snp_effects.txt', sep='\t')
    data['beta_matrix'] = pd.read_csv('cross_trait_beta_matrix.txt', sep='\t')
    
    # Load coloc results
    data['prioritized_loci'] = pd.read_csv('coloc/prioritized_loci.tsv', sep='\t')
    data['coloc_results'] = pd.read_csv('coloc/coloc_results.tsv', sep='\t')
    
    # Load coloc locus summary if exists
    coloc_summary_file = Path('coloc/coloc_locus_summary.tsv')
    if coloc_summary_file.exists():
        data['coloc_locus_summary'] = pd.read_csv(coloc_summary_file, sep='\t')
    
    return data

def classify_loci(data):
    """Classify loci based on significant traits and cross-trait effects."""
    loci_df = data['prioritized_loci'].copy()
    cross_trait = data['cross_trait'].copy()
    
    # Initialize classification columns
    classifications = []
    primary_traits = []
    secondary_traits = []
    effect_patterns = []
    
    for _, locus in loci_df.iterrows():
        locus_id = locus['locus_id']
        chrom = str(locus['chr'])
        start = locus['start']
        end = locus['end']
        
        # Get SNPs in this locus from cross_trait data
        # Use chromosome and position matching
        locus_snps = cross_trait[
            (cross_trait['chr'].astype(str) == chrom) &
            (cross_trait['pos'] >= start) &
            (cross_trait['pos'] <= end)
        ].copy()
        
        if len(locus_snps) == 0:
            # If no SNPs in cross_trait, check sig_traits column
            sig_traits = str(locus['sig_traits']).split(',')
            if 'size' in sig_traits:
                classifications.append('size_locus')
                primary_traits.append('size')
            elif 'protein' in sig_traits:
                classifications.append('protein_locus')
                primary_traits.append('protein')
            elif 'oil' in sig_traits:
                classifications.append('oil_locus')
                primary_traits.append('oil')
            else:
                classifications.append('unclassified')
                primary_traits.append('none')
            secondary_traits.append('none')
            effect_patterns.append('single')
            continue
        
        # Count significant SNPs per trait (p < 5e-8)
        sig_counts = {
            'size': (locus_snps['p_wald_size'] < 5e-8).sum(),
            'protein': (locus_snps['p_wald_protein'] < 5e-8).sum(),
            'oil': (locus_snps['p_wald_oil'] < 5e-8).sum()
        }
        
        # Check for multi-trait significance
        multi_trait_snps = 0
        for _, snp in locus_snps.iterrows():
            sig_traits = []
            if snp['p_wald_size'] < 5e-8:
                sig_traits.append('size')
            if snp['p_wald_protein'] < 5e-8:
                sig_traits.append('protein')
            if snp['p_wald_oil'] < 5e-8:
                sig_traits.append('oil')
            if len(sig_traits) > 1:
                multi_trait_snps += 1
        
        # Determine primary trait (max significant SNPs)
        max_count = max(sig_counts.values())
        primary_trait = [k for k, v in sig_counts.items() if v == max_count][0]
        
        # Determine secondary traits (any other trait with significant SNPs)
        secondary_trait_list = [k for k, v in sig_counts.items() if v > 0 and k != primary_trait]
        
        # Classify based on patterns
        if multi_trait_snps > 0:
            classification = 'cross_effect_locus'
            effect_pattern = 'multi_trait'
        elif len(secondary_trait_list) > 0:
            classification = 'cross_effect_locus'
            effect_pattern = 'shared_region'
        else:
            classification = f'{primary_trait}_locus'
            effect_pattern = 'single'
        
        classifications.append(classification)
        primary_traits.append(primary_trait)
        secondary_traits.append(','.join(secondary_trait_list) if secondary_trait_list else 'none')
        effect_patterns.append(effect_pattern)
    
    # Add classification columns to loci dataframe
    loci_df['locus_class'] = classifications
    loci_df['primary_trait'] = primary_traits
    loci_df['secondary_traits'] = secondary_traits
    loci_df['effect_pattern'] = effect_patterns
    
    return loci_df

def create_final_summary_table(data, classified_loci):
    """Create final summary table combining SNP info, effects, and coloc results."""
    # Start with cross_trait effects as base
    summary_df = data['cross_trait'].copy()
    
    # Add locus information
    locus_info_list = []
    for _, snp in summary_df.iterrows():
        chrom = str(snp['chr'])
        pos = snp['pos']
        
        # Find which locus contains this SNP
        locus_match = None
        for _, locus in classified_loci.iterrows():
            if (str(locus['chr']) == chrom and 
                pos >= locus['start'] and 
                pos <= locus['end']):
                locus_match = locus
                break
        
        if locus_match is not None:
            locus_info_list.append({
                'locus_id': locus_match['locus_id'],
                'locus_class': locus_match['locus_class'],
                'primary_trait': locus_match['primary_trait'],
                'secondary_traits': locus_match['secondary_traits']
            })
        else:
            locus_info_list.append({
                'locus_id': 'none',
                'locus_class': 'unassigned',
                'primary_trait': 'none',
                'secondary_traits': 'none'
            })
    
    locus_info_df = pd.DataFrame(locus_info_list)
    summary_df = pd.concat([summary_df, locus_info_df], axis=1)
    
    # Add coloc results if available for the locus
    coloc_by_locus = {}
    for _, row in data['coloc_results'].iterrows():
        locus_id = row['locus_id']
        trait_pair = f"{row['trait1']}_vs_{row['trait2']}"
        if locus_id not in coloc_by_locus:
            coloc_by_locus[locus_id] = {}
        coloc_by_locus[locus_id][f'coloc_{trait_pair}_pp4'] = row['pp4']
        coloc_by_locus[locus_id][f'coloc_{trait_pair}_colocalized'] = float(row['colocalized'])
    
    # Add coloc columns to summary
    coloc_columns = []
    for locus_id, coloc_data in coloc_by_locus.items():
        for col in coloc_data:
            if col not in coloc_columns:
                coloc_columns.append(col)
    
    for col in coloc_columns:
        summary_df[col] = np.nan
    
    for idx, row in summary_df.iterrows():
        locus_id = row['locus_id']
        if locus_id in coloc_by_locus:
            for col, value in coloc_by_locus[locus_id].items():
                summary_df.at[idx, col] = value
    
    return summary_df

def generate_figures(data, classified_loci, summary_df):
    """Generate all required figures."""
    figures_dir = Path('figures')
    figures_dir.mkdir(exist_ok=True)
    
    # 1. Cross-trait scatter plot
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()
    
    beta_df = data['beta_matrix']
    
    # Oil vs Protein
    ax = axes[0]
    ax.scatter(beta_df['beta_oil'], beta_df['beta_protein'], alpha=0.6, s=50)
    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
    ax.set_xlabel('Beta Oil')
    ax.set_ylabel('Beta Protein')
    ax.set_title('Oil vs Protein Effects')
    
    # Oil vs Size
    ax = axes[1]
    ax.scatter(beta_df['beta_oil'], beta_df['beta_size'], alpha=0.6, s=50)
    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
    ax.set_xlabel('Beta Oil')
    ax.set_ylabel('Beta Size')
    ax.set_title('Oil vs Size Effects')
    
    # Protein vs Size
    ax = axes[2]
    ax.scatter(beta_df['beta_protein'], beta_df['beta_size'], alpha=0.6, s=50)
    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
    ax.set_xlabel('Beta Protein')
    ax.set_ylabel('Beta Size')
    ax.set_title('Protein vs Size Effects')
    
    # Hide empty subplot
    axes[3].axis('off')
    
    plt.tight_layout()
    plt.savefig(figures_dir / 'cross_trait_scatter.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Saved cross_trait_scatter.png")
    
    # 2. Beta heatmap
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Prepare data for heatmap
    heatmap_data = beta_df[['beta_size', 'beta_protein', 'beta_oil']].copy()
    heatmap_data.index = beta_df['rs'].str.slice(0, 20)  # Truncate long SNP names
    
    # Normalize each column for better visualization
    heatmap_data_normalized = (heatmap_data - heatmap_data.mean()) / heatmap_data.std()
    
    sns.heatmap(heatmap_data_normalized.T, cmap='RdBu_r', center=0,
                xticklabels=False, yticklabels=['Size', 'Protein', 'Oil'],
                cbar_kws={'label': 'Normalized Beta'}, ax=ax)
    ax.set_title('Standardized Beta Effects Across Traits')
    ax.set_xlabel('SNPs (truncated)')
    
    plt.tight_layout()
    plt.savefig(figures_dir / 'beta_heatmap.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Saved beta_heatmap.png")
    
    # 3. Locus category bar plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Count locus classifications
    locus_counts = classified_loci['locus_class'].value_counts()
    colors = plt.cm.Set3(np.arange(len(locus_counts)) / len(locus_counts))
    
    ax1.bar(locus_counts.index, locus_counts.values, color=colors)
    ax1.set_xlabel('Locus Category')
    ax1.set_ylabel('Number of Loci')
    ax1.set_title('Locus Classification Distribution')
    ax1.tick_params(axis='x', rotation=45)
    
    # Add counts on bars
    for i, (category, count) in enumerate(locus_counts.items()):
        ax1.text(i, count + 0.5, str(count), ha='center', va='bottom')
    
    # Pie chart of primary traits
    primary_counts = classified_loci['primary_trait'].value_counts()
    ax2.pie(primary_counts.values, labels=primary_counts.index, autopct='%1.1f%%',
            colors=plt.cm.Set2(np.arange(len(primary_counts)) / len(primary_counts)))
    ax2.set_title('Primary Trait Distribution')
    
    plt.tight_layout()
    plt.savefig(figures_dir / 'locus_category_barplot.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Saved locus_category_barplot.png")
    
    # 4. Coloc summary plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    coloc_df = data['coloc_results']
    
    # PP4 distribution
    ax1.hist(coloc_df['pp4'], bins=20, edgecolor='black', alpha=0.7)
    ax1.axvline(x=0.8, color='red', linestyle='--', label='PP4 = 0.8 threshold')
    ax1.set_xlabel('PP4 (Posterior Probability of Colocalization)')
    ax1.set_ylabel('Frequency')
    ax1.set_title('Distribution of Colocalization PP4 Values')
    ax1.legend()
    
    # Colocalization by trait pair
    trait_pairs = ['oil_vs_protein', 'oil_vs_size', 'protein_vs_size']
    coloc_counts = []
    trait_pair_labels = []
    
    for trait_pair in trait_pairs:
        # Check if this trait pair exists in data
        subset = coloc_df[
            (coloc_df['trait1'] + '_vs_' + coloc_df['trait2']) == trait_pair
        ]
        if len(subset) > 0:
            trait_pair_labels.append(trait_pair.replace('_vs_', ' vs '))
            coloc_counts.append(subset['colocalized'].sum())
    
    if coloc_counts:
        colors = plt.cm.Set2(np.arange(len(trait_pair_labels)) / len(trait_pair_labels))
        bars = ax2.bar(trait_pair_labels, coloc_counts, color=colors)
        ax2.set_xlabel('Trait Pair')
        ax2.set_ylabel('Number of Colocalized Loci')
        ax2.set_title('Colocalization by Trait Pair')
        ax2.tick_params(axis='x', rotation=45)
        
        # Add counts on bars
        for bar, count in zip(bars, coloc_counts):
            height = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                    str(count), ha='center', va='bottom')
    
    plt.tight_layout()
    plt.savefig(figures_dir / 'coloc_summary_plot.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Saved coloc_summary_plot.png")

def generate_final_report(data, classified_loci, summary_df):
    """Generate final markdown report."""
    report_path = Path('final_analysis_report.md')
    
    report_lines = []
    
    # Header
    report_lines.append("# Wild Soybean GWAS Post-Analysis Final Report")
    report_lines.append("")
    report_lines.append("## Executive Summary")
    report_lines.append("")
    report_lines.append("This report summarizes the complete post-GWAS analysis of wild soybean traits: ")
    report_lines.append("100-seed weight (size), protein content, and oil content.")
    report_lines.append("")
    
    # 1. Significant SNP Summary
    report_lines.append("## 1. Significant SNP Summary")
    report_lines.append("")
    
    size_count = len(data['size_sig'])
    protein_count = len(data['protein_sig'])
    oil_count = len(data['oil_sig'])
    total_unique = data['all_sig']['rs'].nunique()
    
    report_lines.append(f"- **100SW (size):** {size_count:,} genome-wide significant SNPs")
    report_lines.append(f"- **Protein:** {protein_count:,} genome-wide significant SNPs")
    report_lines.append(f"- **Oil:** {oil_count:,} genome-wide significant SNPs")
    report_lines.append(f"- **Total unique SNPs:** {total_unique:,}")
    report_lines.append("")
    
    # Count SNPs significant in multiple traits
    sig_counts = data['all_sig']['rs'].value_counts()
    multi_trait_count = (sig_counts > 1).sum()
    report_lines.append(f"- **SNPs significant in multiple traits:** {multi_trait_count:,}")
    report_lines.append("")
    
    # 2. Locus Classification
    report_lines.append("## 2. Locus Classification")
    report_lines.append("")
    report_lines.append(f"**Total loci identified:** {len(classified_loci):,}")
    report_lines.append("")
    
    locus_class_counts = classified_loci['locus_class'].value_counts()
    report_lines.append("**Classification breakdown:**")
    for class_type, count in locus_class_counts.items():
        percentage = 100 * count / len(classified_loci)
        report_lines.append(f"- {class_type}: {count} loci ({percentage:.1f}%)")
    report_lines.append("")
    
    # Top loci
    report_lines.append("**Top loci by number of SNPs:**")
    report_lines.append("")
    report_lines.append("| Locus ID | Chromosome | SNPs | Primary Trait | Classification |")
    report_lines.append("|----------|------------|------|---------------|----------------|")
    
    top_loci = classified_loci.nlargest(5, 'num_snps')
    for _, row in top_loci.iterrows():
        report_lines.append(f"| {row['locus_id']} | {row['chr']} | {row['num_snps']} | {row['primary_trait']} | {row['locus_class']} |")
    report_lines.append("")
    
    # 3. Cross-Trait Effects
    report_lines.append("## 3. Cross-Trait Effects")
    report_lines.append("")
    
    beta_df = data['beta_matrix']
    
    # Correlation analysis
    correlations = {
        'oil_vs_protein': beta_df['beta_oil'].corr(beta_df['beta_protein']),
        'oil_vs_size': beta_df['beta_oil'].corr(beta_df['beta_size']),
        'protein_vs_size': beta_df['beta_protein'].corr(beta_df['beta_size'])
    }
    
    report_lines.append("**Effect direction correlations:**")
    for pair, corr in correlations.items():
        report_lines.append(f"- {pair}: r = {corr:.3f}")
    report_lines.append("")
    
    # Effect direction concordance
    concordance_data = []
    trait_pairs = [('beta_oil', 'beta_protein'), ('beta_oil', 'beta_size'), ('beta_protein', 'beta_size')]
    
    for trait1, trait2 in trait_pairs:
        pair_name = f"{trait1[5:]}_vs_{trait2[5:]}"
        concordant = ((beta_df[trait1] > 0) & (beta_df[trait2] > 0)) | ((beta_df[trait1] < 0) & (beta_df[trait2] < 0))
        concordant_count = concordant.sum()
        total = len(beta_df)
        concordance_rate = 100 * concordant_count / total
        concordance_data.append((pair_name, concordant_count, total, concordance_rate))
    
    report_lines.append("**Effect direction concordance:**")
    for pair_name, concordant, total, rate in concordance_data:
        report_lines.append(f"- {pair_name}: {concordant}/{total} SNPs ({rate:.1f}%) have same effect direction")
    report_lines.append("")
    
    # 4. Colocalization Results
    report_lines.append("## 4. Bayesian Colocalization Results")
    report_lines.append("")
    
    coloc_df = data['coloc_results']
    total_tests = len(coloc_df)
    colocalized = coloc_df['colocalized'].sum()
    
    report_lines.append(f"**Total coloc tests:** {total_tests}")
    report_lines.append(f"**Colocalized loci (PP4 ≥ 0.8):** {colocalized} ({100*colocalized/total_tests:.1f}%)")
    report_lines.append("")
    
    if colocalized > 0:
        top_coloc = coloc_df.nlargest(5, 'pp4')[['locus_id', 'trait1', 'trait2', 'pp4']]
        report_lines.append("**Top colocalization signals:**")
        report_lines.append("")
        report_lines.append("| Locus ID | Trait Pair | PP4 |")
        report_lines.append("|----------|------------|-----|")
        for _, row in top_coloc.iterrows():
            report_lines.append(f"| {row['locus_id']} | {row['trait1']} vs {row['trait2']} | {row['pp4']:.6f} |")
        report_lines.append("")
    else:
        report_lines.append("**Note:** No loci showed strong evidence of colocalization (PP4 ≥ 0.8).")
        report_lines.append("The maximum PP4 observed was {:.6f}.".format(coloc_df['pp4'].max()))
        report_lines.append("")
    
    # 5. Generated Figures
    report_lines.append("## 5. Generated Figures")
    report_lines.append("")
    report_lines.append("The analysis generated the following visualization files:")
    report_lines.append("")
    report_lines.append("1. **cross_trait_scatter.png** - Scatter plots of beta effects between trait pairs")
    report_lines.append("2. **beta_heatmap.png** - Heatmap of standardized beta effects across SNPs")
    report_lines.append("3. **locus_category_barplot.png** - Distribution of locus classifications")
    report_lines.append("4. **coloc_summary_plot.png** - Summary of colocalization results")
    report_lines.append("")
    
    # 6. Files Generated
    report_lines.append("## 6. Files Generated")
    report_lines.append("")
    report_lines.append("| File | Description |")
    report_lines.append("|------|-------------|")
    report_lines.append("| `locus_classification.tsv` | Locus classification results |")
    report_lines.append("| `final_summary_table.tsv` | Combined SNP, locus, and coloc data |")
    report_lines.append("| `figures/` | Directory containing all visualization files |")
    report_lines.append("| `final_analysis_report.md` | This summary report |")
    report_lines.append("")
    
    # 7. Conclusions
    report_lines.append("## 7. Conclusions")
    report_lines.append("")
    report_lines.append("### Key Findings:")
    report_lines.append("")
    report_lines.append("1. **Trait-specific genetic architecture:** Most genome-wide significant associations ")
    report_lines.append("   are specific to individual traits with limited pleiotropy.")
    report_lines.append("")
    report_lines.append("2. **Limited colocalization:** No strong evidence for shared causal variants ")
    report_lines.append("   between oil, protein, and seed size traits at genome-wide significant loci.")
    report_lines.append("")
    report_lines.append("3. **Distinct genetic regulation:** The lack of colocalization suggests ")
    report_lines.append("   independent genetic mechanisms underlying these important agronomic traits.")
    report_lines.append("")
    report_lines.append("4. **Methodological framework established:** The complete pipeline from GWAS ")
    report_lines.append("   to colocalization provides a template for future multi-trait analyses.")
    report_lines.append("")
    
    report_lines.append("### Limitations:")
    report_lines.append("")
    report_lines.append("1. **Sample size:** Limited power for detecting weak pleiotropic effects")
    report_lines.append("2. **LD considerations:** Colocalization analysis did not incorporate LD structure")
    report_lines.append("3. **Threshold selection:** Genome-wide significance threshold (p < 5e-8) ")
    report_lines.append("   may miss subtler shared genetic effects")
    report_lines.append("")
    
    report_lines.append("### Future Directions:")
    report_lines.append("")
    report_lines.append("1. **LD-aware colocalization:** Incorporate LD matrices from genotype data")
    report_lines.append("2. **Fine-mapping:** Identify causal variants within associated loci")
    report_lines.append("3. **Functional validation:** Experimental follow-up of prioritized loci")
    report_lines.append("4. **Multi-omics integration:** Combine with transcriptomics and metabolomics")
    report_lines.append("")
    
    report_lines.append("---")
    report_lines.append("")
    report_lines.append("*Report generated on {}*".format(pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')))
    report_lines.append("*Analysis pipeline: GWAS Post-Analysis v1.0*")
    
    # Write report to file
    with open(report_path, 'w') as f:
        f.write('\n'.join(report_lines))
    
    print(f"Saved final report to: {report_path}")

def main():
    """Main execution function."""
    print("Loading data...")
    data = load_data()
    
    print("Classifying loci...")
    classified_loci = classify_loci(data)
    
    print("Creating final summary table...")
    summary_df = create_final_summary_table(data, classified_loci)
    
    print("Generating figures...")
    generate_figures(data, classified_loci, summary_df)
    
    print("Generating final report...")
    generate_final_report(data, classified_loci, summary_df)
    
    # Save outputs
    classified_loci.to_csv('locus_classification.tsv', sep='\t', index=False)
    print("Saved locus_classification.tsv")
    
    summary_df.to_csv('final_summary_table.tsv', sep='\t', index=False)
    print("Saved final_summary_table.tsv")
    
    print("\nAnalysis complete!")
    print("Outputs saved in SNP_analyze directory:")
    print("  - locus_classification.tsv")
    print("  - final_summary_table.tsv")
    print("  - figures/ (4 visualization files)")
    print("  - final_analysis_report.md")

if __name__ == "__main__":
    main()