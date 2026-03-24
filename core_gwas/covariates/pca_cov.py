# Generalized from soybean project-specific path layout.
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re

def normalize_id(iid):
    """
    将 GDW094 转换为 Wild94, GDC001 转换为 Cultivar1
    """
    match = re.match(r'([a-zA-Z]+)([0-9]+)', str(iid))
    if not match: return iid
    prefix, num = match.groups()
    num_int = int(num) # 自动移除前导 0
    mapping = {'GDC': 'Cultivar', 'GDL': 'Landrace', 'GDW': 'Wild'}
    new_prefix = mapping.get(prefix, prefix)
    return f"{new_prefix}{num_int}"

def load_data():
    # 1. 读取 PCA 结果
    file_prefix = 'soybean_pca'
    eigenval = pd.read_csv(f'{file_prefix}.eigenval', header=None, names=['val'])
    
    df_pca = pd.read_csv(f'{file_prefix}.eigenvec', sep=r'\s+')
    df_pca.columns = [c.replace('#', '') for c in df_pca.columns]
    
    # 2. 读取表型数据 (TSV格式)
    pheno_path = 'data/input/workflow/phenotype_long.tsv'
    df_pheno = pd.read_csv(pheno_path, sep='\t')
    
    # 3. ID 归一化与合并
    # 假设表型文件的 ID 列名为 'ID' 或 'Accession'，请根据实际情况微调
    # 这里我们统一创建一个 Match_ID
    df_pca['Match_ID'] = df_pca['IID'].apply(normalize_id)
    # 假设 phenotype_long.tsv 里的 ID 列叫 'ID'
    df_pheno['Match_ID'] = df_pheno['ID'].astype(str).str.strip() 
    
    combined = pd.merge(df_pca, df_pheno, on='Match_ID', how='inner')
    
    # 4. 自动分类标签
    def assign_group(iid):
        if 'GDW' in iid or 'Wild' in iid: return 'Wild (野生)'
        if 'GDL' in iid or 'Landrace' in iid: return 'Landrace (地方种)'
        if 'GDC' in iid or 'Cultivar' in iid: return 'Cultivar (栽培)'
        return 'Other'
    
    combined['Group'] = combined['IID'].apply(assign_group)
    return combined, eigenval

def plot_enhanced_pca(df, eigenval):
    total_ev = eigenval['val'].sum()
    ev_per = (eigenval['val'] / total_ev * 100).round(2)

    plt.figure(figsize=(10, 8))
    sns.set_style("whitegrid")
    
    # 定义进化色系
    colors = {'Wild (野生)': '#d62728', 'Landrace (地方种)': '#ff7f0e', 'Cultivar (栽培)': '#1f77b4', 'Other': '#7f7f7f'}
    
    # 绘图
    scatter = sns.scatterplot(data=df, x='PC1', y='PC2', hue='Group', 
                               palette=colors, s=80, alpha=0.8, edgecolor='w')
    
    plt.xlabel(f'PC1 ({ev_per[0]}%)', fontsize=12)
    plt.ylabel(f'PC2 ({ev_per[1]}%)', fontsize=12)
    plt.title('Soybean Population Structure: Domestication Trajectory', fontsize=15)
    
    # 添加样本量标注
    info_text = f"Matched Samples: {len(df)}"
    plt.text(0.95, 0.05, info_text, transform=plt.gca().transAxes, 
             fontsize=10, horizontalalignment='right',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

    plt.savefig('soybean_evolution_pca.png', dpi=300, bbox_inches='tight')
    print(f"✓ PCA 演化分布图已生成: soybean_evolution_pca.png")

def plot_trait_comparison(df):
    # 动态检测表型列
    traits = ['Protein', 'Oil', '100SW']
    # 过滤掉不存在的列
    available_traits = [t for t in traits if t in df.columns]
    
    if not available_traits:
        print("未在合并数据中发现 Protein/Oil/100SW 列，请检查 TSV 列名")
        return

    plt.figure(figsize=(5 * len(available_traits), 6))
    colors = {'Wild (野生)': '#d62728', 'Landrace (地方种)': '#ff7f0e', 'Cultivar (栽培)': '#1f77b4', 'Other': '#7f7f7f'}

    for i, trait in enumerate(available_traits):
        plt.subplot(1, len(available_traits), i + 1)
        sns.boxplot(x='Group', y=trait, data=df, palette=colors)
        sns.stripplot(x='Group', y=trait, data=df, color='black', size=2, alpha=0.3)
        plt.title(f'{trait} Distribution')
    
    plt.tight_layout()
    plt.savefig('soybean_traits_by_group.png', dpi=300)
    print(f"✓ 表型亚群对比图已生成: soybean_traits_by_group.png")

if __name__ == "__main__":
    combined_df, eigenval_df = load_data()
    print(f"合并后样本总数: {len(combined_df)}")
    if len(combined_df) > 0:
        plot_enhanced_pca(combined_df, eigenval_df)
        plot_trait_comparison(combined_df)
    else:
        print("❌ 错误：合并样本数为 0，请检查 ID 匹配逻辑。")