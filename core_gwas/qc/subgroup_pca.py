# Generalized from soybean project-specific path layout.
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re
from scipy import stats


def normalize_id(iid):
    """将 GDW094 转换为 Wild94, GDC001 转换为 Cultivar1"""
    match = re.match(r'([a-zA-Z]+)([0-9]+)', str(iid).strip())
    if not match:
        return str(iid).strip()
    prefix, num = match.groups()
    num_int = int(num)  # 自动移除前导0
    mapping = {'GDC': 'Cultivar', 'GDL': 'Landrace', 'GDW': 'Wild'}
    new_prefix = mapping.get(prefix, prefix)
    return f"{new_prefix}{num_int}"


def load_data(pca_prefix, pheno_path, id_col, trait_cols):
    """
    针对 PLINK2 格式优化的数据加载函数
    """
    print(f"正在加载 PCA 文件: {pca_prefix}.eigenvec")
    
    # 1. 读取 PCA 文件 (处理可能存在的 #IID 标题)
    # sep=r'\s+' 兼容空格或制表符
    try:
        df_pca_raw = pd.read_csv(f'{pca_prefix}.eigenvec', sep=r'\s+', comment='#', header=None)
    except Exception as e:
        print(f"❌ 读取 PCA 文件失败: {e}")
        return None, None

    # --- 自动对位 ID 列 ---
    # 观察数据发现第一列(index 0)是 GDC001 这种 ID
    # 如果第一列包含数字（且不含字母），则可能是数值列；如果含字母如 'GDC'，则是 ID
    first_cell = str(df_pca_raw.iloc[0, 0])
    if any(c.isalpha() for c in first_cell):
        id_idx = 0
    else:
        id_idx = 1 # 处理 FID IID 都在的情况

    df_pca = pd.DataFrame()
    df_pca['Match_ID'] = df_pca_raw.iloc[:, id_idx].astype(str).str.strip().str.upper()
    df_pca['PC1'] = pd.to_numeric(df_pca_raw.iloc[:, id_idx + 1], errors='coerce')
    df_pca['PC2'] = pd.to_numeric(df_pca_raw.iloc[:, id_idx + 2], errors='coerce')
    df_pca['PC3'] = pd.to_numeric(df_pca_raw.iloc[:, id_idx + 3], errors='coerce')
    df_pca['PC4'] = pd.to_numeric(df_pca_raw.iloc[:, id_idx + 4], errors='coerce')

    print(f"✓ PCA 加载成功，Match_ID 示例: {df_pca['Match_ID'].iloc[0]}")

# --- 重点修改表型文件读取 ---
    print(f"正在加载表型文件: {pheno_path}")
    try:
        # 显式指定制表符分隔，并去掉列名的前后空格
        df_pheno = pd.read_csv(pheno_path, sep='\t')
        df_pheno.columns = df_pheno.columns.str.strip() 
        
        # 调试：打印出程序实际看到的列名，确保 'ID' 在里面
        print(f"✅ 成功读取表型列: {df_pheno.columns.tolist()}")
    except Exception as e:
        print(f"❌ 读取失败: {e}")
        return None, None

    if id_col not in df_pheno.columns:
        print(f"❌ 错误：在文件中找不到列名 '{id_col}'")
        return None, None

    # 处理匹配 ID
    df_pheno['Match_ID'] = df_pheno[id_col].astype(str).str.strip().str.upper()
    
    # 3. 合并
    combined_df = pd.merge(df_pca, df_pheno, on='Match_ID', how='inner')
    
    # 4. 读取特征值 (eigenval)
    try:
        eigenval_df = pd.read_csv(f'{pca_prefix}.eigenval', header=None)
    except:
        eigenval_df = pd.DataFrame()

    return combined_df, eigenval_df

    
    # 4. 亚群分类
    def assign_group(iid):
        iid_str = str(iid).upper()
        if 'GDW' in iid_str: return 'Wild (野生)'
        if 'GDL' in iid_str: return 'Landrace (地方种)'
        if 'GDC' in iid_str: return 'Cultivar (栽培)'
        return 'Other'

    combined['Group'] = combined['IID_orig'].apply(assign_group)

    # 5. 表型QC
    trait_list = [col for col in trait_cols.split(',') if col.strip() in combined.columns]
    if trait_list:
        combined['OutlierFlag'] = False
        for trait in trait_list:
            mean_val = combined[trait].mean()
            std_val = combined[trait].std()
            threshold = 3.0
            outliers = (np.abs(combined[trait] - mean_val) > threshold * std_val)
            combined.loc[outliers, 'OutlierFlag'] = True
            print(f" - {trait}: {outliers.sum()} 个异常值 (|Z| > {threshold}σ)")
    
    return combined, eigenval


def plot_enhanced_pca(df, eigenval):
    """绘制进化轨迹PCA图"""
    total_ev = eigenval['val'].sum()
    ev_per = (eigenval['val'] / total_ev * 100).round(2)
    
    plt.figure(figsize=(12, 8))
    sns.set_style("whitegrid")
    
    colors = {
        'Wild (野生)': '#d62728', 
        'Landrace (地方种)': '#ff7f0e', 
        'Cultivar (栽培)': '#1f77b4', 
        'Other': '#7f7f7f'
    }
    
    sns.scatterplot(
        data=df, x='PC1', y='PC2', hue='Group', 
        palette=colors, s=100, alpha=0.7, edgecolor='black', linewidth=0.5
    )
    
    plt.xlabel(f'PC1 ({ev_per[0]}%)', fontsize=13)
    plt.ylabel(f'PC2 ({ev_per[1]}%)', fontsize=13)
    plt.title('大豆群体结构: 驯化轨迹\n(三个亚群进化关系)', fontsize=14, fontweight='bold')

    plt.tight_layout()
    plt.savefig('soybean_evolution_pca.png', dpi=300, bbox_inches='tight')
    print("✓ PCA进化图已生成: soybean_evolution_pca.png")
    plt.close()


def plot_trait_comparison(df, trait_cols):
    """绘制表型亚群对比"""
    traits = [t.strip() for t in trait_cols.split(',')]
    available_traits = [t for t in traits if t in df.columns]
    
    if not available_traits:
        print("⚠️ 未找到指定的表型列")
        return

    fig, axes = plt.subplots(1, len(available_traits), figsize=(5 * len(available_traits), 6))
    if len(available_traits) == 1:
        axes = [axes]
    
    colors = {
        'Wild (野生)': '#d62728', 
        'Landrace (地方种)': '#ff7f0e', 
        'Cultivar (栽培)': '#1f77b4', 
        'Other': '#7f7f7f'
    }
    
    for i, trait in enumerate(available_traits):
        ax = axes[i]
        sns.boxplot(x='Group', y=trait, data=df, palette=colors, ax=ax)
        sns.stripplot(x='Group', y=trait, data=df, color='black', size=3, alpha=0.4, ax=ax)
        ax.set_title(f'{trait} 亚群分布', fontweight='bold')
        
        # ANOVA 检验
        groups = [g for g in df['Group'].unique() if g != 'Other']
        group_data = [df[df['Group'] == g][trait].dropna().values for g in groups]
        if len(group_data) > 1:
            _, p_val = stats.f_oneway(*group_data)
            ax.text(0.5, 0.95, f'ANOVA p={p_val:.2e}', transform=ax.transAxes, 
                    ha='center', va='top', bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.3))

    plt.tight_layout()
    plt.savefig('soybean_traits_by_group.png', dpi=300)
    print("✓ 表型对比图已生成: soybean_traits_by_group.png")
    plt.close()


if __name__ == "__main__":
    pca_prefix = 'soybean_pca'
    pheno_path = 'data/input/workflow/phenotype_long.tsv'
    id_col = 'ID'  # 保持为 ID
    trait_cols = ['Protein', 'Oil', '100SW']

    combined_df, eigenval_df = load_data(pca_prefix, pheno_path, id_col, trait_cols)
    
    if combined_df is not None and len(combined_df) > 0:
        plot_enhanced_pca(combined_df, eigenval_df)
        plot_trait_comparison(combined_df, trait_cols)
        print("\n✅ 分析完成！")
    else:
        print("\n❌ 分析失败，请检查输入文件和参数")