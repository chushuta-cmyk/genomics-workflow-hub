import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats

def load_phenotype_data(pheno_file, id_col='ID'):
    """Loads phenotype data and assigns groups based on ID prefixes"""
    print(f"Loading phenotype file: {pheno_file}")
    try:
        # 1. Load data
        df = pd.read_csv(pheno_file, sep='\t')
        df.columns = df.columns.str.strip()
        
        # 2. Calculate mean across years to remove year effect
        traits = ['Protein', 'Oil', '100SW']
        # Filter traits that actually exist in the file
        existing_traits = [t for t in traits if t in df.columns]
        
        pheno_avg = df.groupby(id_col)[existing_traits].mean().reset_index()
        
        # 3. Assign Sub-populations based on ID prefix
        def assign_group(id_val):
            val = str(id_val).upper()
            if 'GDW' in val: return 'Wild'
            elif 'GDL' in val: return 'Landrace'
            elif 'GDC' in val: return 'Cultivar'
            return 'Other'
        
        pheno_avg['Group'] = pheno_avg[id_col].apply(assign_group)
        
        print(f"✓ Phenotype data processed. Total unique IDs: {len(pheno_avg)}")
        print("Group distribution:")
        print(pheno_avg['Group'].value_counts())
        
        return pheno_avg, existing_traits
    except Exception as e:
        print(f"❌ Error loading phenotype file: {e}")
        return None, []

def plot_correlation_matrix(df, traits):
    """Generates the PairGrid correlation plot"""
    print("Generating correlation grid...")
    sns.set_style("white")
    
    def corrfunc(x, y, **kws):
        mask = ~np.isnan(x) & ~np.isnan(y)
        if sum(mask) < 2: return
        r, p = stats.pearsonr(x[mask], y[mask])
        ax = plt.gca()
        ax.annotate(f"r = {r:.2f}\np = {p:.1e}",
                    xy=(0.1, 0.9), xycoords=ax.transAxes,
                    fontsize=10, fontweight='bold', color='red')

    g = sns.PairGrid(df[traits + ['Group']].dropna(subset=traits), aspect=1.2)
    g.map_upper(corrfunc)
    g.map_lower(sns.regplot, 
                scatter_kws={'alpha':0.4, 's':20, 'color':'teal'},
                line_kws={'color':'red', 'lw':1.5})
    g.map_diag(sns.histplot, kde=True, color='steelblue', edgecolor='w')
    
    g.fig.suptitle('Soybean Traits Relationship Analysis', y=1.02, fontsize=14)
    plt.savefig('soybean_traits_correlation.png', dpi=300, bbox_inches='tight')
    plt.close()

def plot_subgroup_boxplots(df, traits):
    """Generates boxplots with ANOVA stats for sub-populations"""
    print("Generating subgroup boxplots...")
    n_traits = len(traits)
    fig, axes = plt.subplots(1, n_traits, figsize=(5 * n_traits, 6))
    if n_traits == 1: axes = [axes]
    
    colors = {'Wild': '#d62728', 'Landrace': '#ff7f0e', 'Cultivar': '#1f77b4'}
    order = ['Wild', 'Landrace', 'Cultivar']
    # Only plot groups that actually exist in the data
    plot_order = [g for g in order if g in df['Group'].unique()]

    for i, trait in enumerate(traits):
        ax = axes[i]
        sns.boxplot(data=df, x='Group', y=trait, palette=colors, ax=ax, order=plot_order, showfliers=False)
        sns.stripplot(data=df, x='Group', y=trait, color='black', alpha=0.3, ax=ax, order=plot_order)
        
        # ANOVA
        group_data = [df[df['Group'] == g][trait].dropna() for g in plot_order]
        if len(group_data) > 1:
            _, p_val = stats.f_oneway(*group_data)
            ax.set_title(f'{trait}\n(ANOVA p={p_val:.2e})', fontweight='bold')
        else:
            ax.set_title(f'{trait}', fontweight='bold')
        
        ax.set_xlabel('')

    plt.tight_layout()
    plt.savefig('soybean_traits_boxplot_subgroups.png', dpi=300)
    plt.close()

if __name__ == "__main__":
    # Settings
    input_file = 'phenotype_long.tsv'
    id_column = 'ID'
    
    # 1. Process Data
    pheno_avg, trait_list = load_phenotype_data(input_file, id_column)
    
    if pheno_avg is not None:
        # 2. Run Visualizations
        plot_correlation_matrix(pheno_avg, trait_list)
        plot_subgroup_boxplots(pheno_avg, trait_list)
        
        # 3. Protein vs Oil Scatter (Categorized by Group)
        if 'Protein' in trait_list and 'Oil' in trait_list:
            plt.figure(figsize=(8, 6))
            sns.scatterplot(data=pheno_avg, x='Protein', y='Oil', hue='Group', 
                            palette={'Wild': '#d62728', 'Landrace': '#ff7f0e', 'Cultivar': '#1f77b4'}, 
                            s=100, alpha=0.7)
            plt.title('Protein vs Oil Trade-off by Sub-population')
            plt.grid(True, linestyle='--', alpha=0.5)
            plt.savefig('protein_vs_oil_subgroup_scatter.png', dpi=300)
            plt.close()

        # 4. Save Summary Statistics
        summary = pheno_avg.groupby('Group')[trait_list].agg(['mean', 'std', 'count'])
        summary.to_csv('phenotype_summary_by_group.csv')
        
        print("\n" + "="*40)
        print("✅ Phenotype Analysis Complete!")
        print("Generated: soybean_traits_correlation.png")
        print("Generated: soybean_traits_boxplot_subgroups.png")
        print("Generated: protein_vs_oil_subgroup_scatter.png")
        print("Generated: phenotype_summary_by_group.csv")
        print("="*40)