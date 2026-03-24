# Generalized from soybean project-specific path layout.
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os
import re

def clean_id(val):
    if pd.isna(val): return ""
    s = str(val).strip()
    match = re.match(r'([a-zA-Z]+)([0-9]+)', s)
    if match:
        prefix, num = match.groups()
        return f"{prefix.upper()}{int(num)}"
    return s.upper()

def load_and_merge():
    # 1. Load PCA
    pca_file = 'pca_data.eigenvec'
    df_pca_raw = pd.read_csv(pca_file, sep='\s+')
    df_pca_raw.columns = [c.replace('#', '') for c in df_pca_raw.columns]
    df_pca = df_pca_raw[['IID', 'PC1', 'PC2', 'PC3']].copy()
    df_pca.columns = ['ID_pca', 'PC1', 'PC2', 'PC3']
    
    # 2. Load Phenotype
    pheno_path = 'data/input/data/phenotype_original.tsv'
    df_pheno = pd.read_csv(pheno_path, sep='\t')
    
    # 3. Normalize IDs
    df_pca['match_id'] = df_pca['ID_pca'].apply(clean_id)
    id_col_name = df_pheno.columns[1] 
    df_pheno['match_id'] = df_pheno[id_col_name].apply(clean_id)
    
    # 4. Merge
    combined = pd.merge(df_pca, df_pheno, on='match_id', how='inner')
    
    # 5. Assign Groups (English labels)
    def assign_group(mid):
        if 'GDC' in mid: return 'Cultivar'
        if 'GDL' in mid: return 'Landrace'
        return 'Other'
    
    combined['Group'] = combined['match_id'].apply(assign_group)
    combined = combined[combined['Group'] != 'Other']
    
    return combined

def plot_subgroups(df):
    plt.figure(figsize=(10, 8))
    sns.set_style("whitegrid")
    colors = {'Cultivar': '#1f77b4', 'Landrace': '#ff7f0e'}
    
    sns.scatterplot(data=df, x='PC1', y='PC2', hue='Group', 
                    palette=colors, s=100, alpha=0.8, edgecolor='w')
    
    plt.title('Soybean Subgroup Structure (Cultivar vs Landrace)', fontsize=15)
    plt.xlabel('Principal Component 1')
    plt.ylabel('Principal Component 2')
    plt.savefig('subgroup_pca_plot.png', dpi=300)
    print("✓ PCA plot saved: subgroup_pca_plot.png")

def plot_traits(df):
    traits = ['Protein', 'Oil', '100SW']
    available_traits = [t for t in traits if t in df.columns]
    if not available_traits: return

    plt.figure(figsize=(15, 6))
    colors = {'Cultivar': '#1f77b4', 'Landrace': '#ff7f0e'}

    for i, trait in enumerate(available_traits):
        plt.subplot(1, 3, i + 1)
        sns.boxplot(x='Group', y=trait, data=df, palette=colors)
        sns.stripplot(x='Group', y=trait, data=df, color='black', size=2, alpha=0.3)
        
        # T-test for significance
        g1 = df[df['Group'] == 'Cultivar'][trait].dropna()
        g2 = df[df['Group'] == 'Landrace'][trait].dropna()
        t_stat, p_val = stats.ttest_ind(g1, g2)
        plt.title(f'{trait}\np-val: {p_val:.2e}')
    
    plt.tight_layout()
    plt.savefig('subgroup_trait_comparison.png', dpi=300)
    print("✓ Trait comparison plot saved: subgroup_trait_comparison.png")

if __name__ == "__main__":
    try:
        combined_data = load_and_merge()
        if len(combined_data) > 0:
            # FUNCTIONS UNCOMMENTED HERE:
            plot_subgroups(combined_data)
            plot_traits(combined_data)
            print(f"✓ Success! Processed {len(combined_data)} records.")
        else:
            print("❌ Error: No samples matched.")
    except Exception as e:
        print(f"❌ Execution Error: {e}")