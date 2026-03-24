# Generalized from soybean project-specific path layout.
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

def load_and_merge_data(pca_file, pheno_file, id_col_name='ID', trait_cols=['Protein', 'Oil', '100SW']):
    """Fixed data loading and preprocessing function"""
    print("="*60)
    print("STEP 1: Reading PCA file")
    print("="*60)
    try:
        # Load PCA eigenvectors (.eigenvec)
        df_pca = pd.read_csv(pca_file, sep=r'\s+', comment='#', header=None)
        df_pca_clean = pd.DataFrame()
        
        # In PLINK, typically col 0 is FID and col 1 is IID. 
        # Using index 0 based on your previous successful Match_ID: GDC001
        df_pca_clean['IID_orig'] = df_pca.iloc[:, 0].astype(str).str.strip()
        df_pca_clean['Match_ID'] = df_pca_clean['IID_orig'].str.upper()
        
        # Extract first 10 PCs
        for i in range(10):
            df_pca_clean[f'PC{i+1}'] = pd.to_numeric(df_pca.iloc[:, i+1], errors='coerce')
        print(f"✓ PCA loaded successfully. Samples: {len(df_pca_clean)}")
    except Exception as e:
        print(f"❌ PCA loading failed: {e}")
        return None, None

    print("\n" + "="*60)
    print("STEP 2: Reading Phenotype file")
    print("="*60)
    try:
        # Using sep='\t' for your .tsv file
        df_pheno = pd.read_csv(pheno_file, sep='\t')
        df_pheno.columns = df_pheno.columns.str.strip()
        
        if id_col_name not in df_pheno.columns:
            print(f"❌ ID column '{id_col_name}' not found. Available: {df_pheno.columns.tolist()}")
            return None, None
            
        df_pheno['Match_ID'] = df_pheno[id_col_name].astype(str).str.strip().str.upper()
        print(f"✓ Phenotype loaded successfully. Samples: {len(df_pheno)}")
    except Exception as e:
        print(f"❌ Phenotype loading failed: {e}")
        return None, None

    print("\n" + "="*60)
    print("STEP 3: Merging Data")
    print("="*60)
    combined = pd.merge(df_pca_clean, df_pheno, on='Match_ID', how='inner')
    print(f"✓ Merge successful: {len(combined)} overlapping samples found.")

    if len(combined) == 0:
        print("❌ Warning: No matching IDs found between PCA and Phenotype files.")
        return None, None

    print("\n" + "="*60)
    print("STEP 4: Sub-population Classification")
    print("="*60)
    def assign_group(iid_str):
        iid_upper = str(iid_str).upper()
        if 'GDW' in iid_upper: return 'Wild'
        elif 'GDL' in iid_upper: return 'Landrace'
        elif 'GDC' in iid_upper: return 'Cultivar'
        else: return 'Other'
    
    combined['Group'] = combined['IID_orig'].apply(assign_group)
    
    # Load eigenvalues for Variance Explained calculations
    try:
        eigenval_file = pca_file.replace('.eigenvec', '.eigenval')
        # Fixing the KeyError: 'val' by providing column names
        eigenval_df = pd.read_csv(eigenval_file, header=None, names=['val'])
    except:
        print("⚠️ .eigenval file not found. PC labels will not show percentages.")
        eigenval_df = pd.DataFrame()

    return combined, eigenval_df

def plot_pca_by_subpop(combined_df):
    """Generates a 3-panel PCA plot split by group"""
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    colors = {'Wild': '#d62728', 'Landrace': '#ff7f0e', 'Cultivar': '#1f77b4'}
    groups = ['Wild', 'Landrace', 'Cultivar']
    
    for idx, group in enumerate(groups):
        df_sub = combined_df[combined_df['Group'] == group]
        ax = axes[idx]
        ax.scatter(df_sub['PC1'], df_sub['PC2'], c=colors.get(group, 'gray'), s=60, alpha=0.6, edgecolor='w')
        ax.set_title(f'{group} (N={len(df_sub)})')
        ax.set_xlabel('PC1')
        ax.set_ylabel('PC2')
        ax.grid(True, linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    plt.savefig('soybean_pca_by_subpop.png', dpi=300)
    plt.close()

def plot_combined_pca(combined_df, eigenval_df):
    """Generates a single combined PCA plot with legend"""
    plt.figure(figsize=(10, 7))
    colors = {'Wild': '#d62728', 'Landrace': '#ff7f0e', 'Cultivar': '#1f77b4'}
    
    sns.scatterplot(data=combined_df, x='PC1', y='PC2', hue='Group', palette=colors, s=80, alpha=0.8)
    
    if not eigenval_df.empty:
        total_ev = eigenval_df['val'].sum()
        pc1_var = (eigenval_df['val'].iloc[0] / total_ev) * 100
        pc2_var = (eigenval_df['val'].iloc[1] / total_ev) * 100
        plt.xlabel(f'PC1 ({pc1_var:.2f}%)')
        plt.ylabel(f'PC2 ({pc2_var:.2f}%)')
    else:
        plt.xlabel('PC1')
        plt.ylabel('PC2')
    
    plt.title('Soybean Population Structure (PCA)')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, linestyle='--', alpha=0.3)
    plt.tight_layout()
    plt.savefig('soybean_pca_combined.png', dpi=300)
    plt.close()

def plot_trait_comparison(combined_df, trait_cols):
    """Generates boxplots for phenotype traits across groups"""
    available_traits = [t for t in trait_cols if t in combined_df.columns]
    if not available_traits: return
    
    fig, axes = plt.subplots(1, len(available_traits), figsize=(6*len(available_traits), 6))
    if len(available_traits) == 1: axes = [axes]
    
    colors = {'Wild': '#d62728', 'Landrace': '#ff7f0e', 'Cultivar': '#1f77b4'}
    order = ['Wild', 'Landrace', 'Cultivar']
    
    for i, trait in enumerate(available_traits):
        ax = axes[i]
        sns.boxplot(data=combined_df, x='Group', y=trait, palette=colors, ax=ax, order=order)
        sns.stripplot(data=combined_df, x='Group', y=trait, color='black', size=3, alpha=0.3, ax=ax, order=order)
        ax.set_title(f'{trait} across Sub-populations')
        
    plt.tight_layout()
    plt.savefig('soybean_traits_by_subpop.png', dpi=300)
    plt.close()

if __name__ == "__main__":
    # Configuration
    pca_file = 'soybean_pca.eigenvec'
    pheno_file = 'data/input/workflow/phenotype_long.tsv'
    id_col_name = 'ID'
    trait_cols = ['Protein', 'Oil', '100SW']

    # Execution
    combined_df, eigenval_df = load_and_merge_data(pca_file, pheno_file, id_col_name, trait_cols)
    
    if combined_df is not None:
        plot_pca_by_subpop(combined_df)
        plot_combined_pca(combined_df, eigenval_df)
        plot_trait_comparison(combined_df, trait_cols)
        print("\n✅ Analysis Complete. Files generated:")
        print(" - soybean_pca_by_subpop.png")
        print(" - soybean_pca_combined.png")
        print(" - soybean_traits_by_subpop.png")