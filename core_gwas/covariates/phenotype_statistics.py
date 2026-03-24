import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, gaussian_kde
from pingouin import partial_corr
# 1. 核心逻辑：合并表型与样本信息 (解决 KeyError: 'Group')
def prepare_data(pheno_file, sample_info_file):
    # 读取表型数据
    df_pheno = pd.read_csv(pheno_file, sep='\t')
    
    # 读取样本信息 (假设该文件包含 'ID' 和 'Group' 列)
    # 这里的 sample_info_file 通常是你做 PCA 或 GWAS 时的样本列表
    df_group = pd.read_csv(sample_info_file, sep='\t')
    
    # 合并数据
    # 使用 inner join 确保只有既有表型又有分组信息的样本被保留
    combined_df = pd.merge(df_pheno, df_group[['ID', 'Group']], on='ID')
    return combined_df

# 2. 统计逻辑：偏相关分析
def run_statistics(df):
    if all(col in df.columns for col in ['Protein', 'Oil', '100SW']):
        res = partial_corr(data=df, x='Protein', y='Oil', covar='100SW')
        print("\n--- 偏相关分析报告 (排除 Size 效应) ---")
        print(res)

def plot_protein_oil_scatter_optimized(combined_df):
    """
    优化后的散点图：展示 Protein vs Oil，
    点大小代表 100SW，并计算各亚群特异性相关系数。
    """
    print("正在生成优化版 Protein vs Oil 散点图...")
    fig, ax = plt.subplots(figsize=(14, 10))
    
    #
    colors = {
        'Wild (野生)': '#d62728',
        'Landrace (地方种)': '#ff7f0e', 
        'Cultivar (栽培)': '#1f77b4'
    }
    
    stats_text = []
    
    for group, color in colors.items():
        df_sub = combined_df[combined_df['Group'] == group].dropna(subset=['Protein', 'Oil', '100SW'])
        if df_sub.empty: continue
        
        # 归一化点大小 (30 到 300 范围)
        size_min, size_max = combined_df['100SW'].min(), combined_df['100SW'].max()
        sizes = (df_sub['100SW'] - size_min) / (size_max - size_min) * 270 + 30
        
        ax.scatter(
            df_sub['Protein'], df_sub['Oil'],
            c=color, s=sizes, alpha=0.5,
            edgecolor='white', linewidth=0.8, label=f"{group} (N={len(df_sub)})"
        )
        
        # 计算亚群内的统计相关性
        r_po, _ = pearsonr(df_sub['Protein'], df_sub['Oil'])
        r_ps = df_sub['Protein'].corr(df_sub['100SW'])
        r_os = df_sub['Oil'].corr(df_sub['100SW'])
        
        stats_text.append(f"• {group}:\n  r(P-O)={r_po:.2f}\n  r(P-Size)={r_ps:.2f}\n  r(O-Size)={r_os:.2f}")

    # 绘制整体趋势线
    sns.regplot(data=combined_df, x='Protein', y='Oil', scatter=False, 
                color='black', order=2, ax=ax, line_kws={"ls":"--", "lw":2, "label":"整体趋势"})
    
    # 在图表侧面添加统计文字框
    plt.text(1.02, 0.5, "\n\n".join(stats_text), transform=ax.transAxes, 
             fontsize=11, verticalalignment='center', bbox=dict(facecolor='white', alpha=0.8))
    
    ax.set_xlabel('Protein Content (%)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Oil Content (%)', fontsize=12, fontweight='bold')
    ax.set_title('Soybean Protein-Oil Trade-off modulated by Seed Weight (100SW)', fontsize=15)
    ax.legend(loc='upper right', frameon=True)
    
    plt.tight_layout()
    plt.savefig('soybean_optimized_scatter.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("✓ 已保存: soybean_optimized_scatter.png")

def plot_trait_distribution_subpop(combined_df, trait_cols):
    """绘制亚群分布密度图"""
    print("正在生成分布密度图...")
    available_traits = [t for t in trait_cols if t in combined_df.columns]
    fig, axes = plt.subplots(1, len(available_traits), figsize=(5*len(available_traits), 5))
    
    groups = ['Wild (野生)', 'Landrace (地方种)', 'Cultivar (栽培)']
    colors = ['#d62728', '#ff7f0e', '#1f77b4']
    
    for i, trait in enumerate(available_traits):
        ax = axes[i] if len(available_traits) > 1 else axes
        for g_idx, group in enumerate(groups):
            data = combined_df[combined_df['Group'] == group][trait].dropna()
            if not data.empty:
                sns.kdeplot(data, ax=ax, label=group, color=colors[g_idx], fill=True, alpha=0.3)
        
        ax.set_title(f'{trait} Distribution', fontweight='bold')
        ax.legend()

    plt.tight_layout()
    plt.savefig('soybean_trait_distribution_by_subpop.png', dpi=300)
    plt.close()

def plot_correlation_heatmap_subpop(combined_df, trait_cols):
    """亚群相关性热图"""
    print("正在生成相关性热图...")
    groups = ['Wild (野生)', 'Landrace (地方种)', 'Cultivar (栽培)']
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    for idx, group in enumerate(groups):
        df_sub = combined_df[combined_df['Group'] == group][trait_cols].dropna()
        if len(df_sub) > 1:
            corr = df_sub.corr()
            sns.heatmap(corr, annot=True, cmap='RdBu_r', center=0, ax=axes[idx], fmt=".2f")
            axes[idx].set_title(f'{group} Correlation')
    
    plt.tight_layout()
    plt.savefig('soybean_correlation_heatmap_by_subpop.png', dpi=300)
    plt.close()

# --- 执行入口 ---
if __name__ == "__main__":
    # 请确保文件名与你的 tsv 一致
    analyze_phenotypes('phenotype_long.tsv')