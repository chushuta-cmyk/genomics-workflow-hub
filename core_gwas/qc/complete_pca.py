import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def load_pca_data(prefix):
    """
    健壮地加载 PLINK2 生成的 PCA 结果
    """
    # 1. 读取特征值 (Eigenvalues)
    eigenval = pd.read_csv(f'{prefix}.eigenval', header=None, names=['val'])
    
    # 2. 读取特征向量 (Eigenvectors)
    # PLINK2 默认带标题行且第一列是 #FID
    try:
        df = pd.read_csv(f'{prefix}.eigenvec', sep=r'\s+')
        # 移除 # 符号并统一列名
        df.columns = [c.replace('#', '') for c in df.columns]
    except Exception as e:
        print(f"读取失败，尝试原始模式: {e}")
        df = pd.read_csv(f'{prefix}.eigenvec', sep=r'\s+', header=None)
    
    # 3. 强制转换 PC 列为数值型，处理可能存在的字符粘连
    pc_cols = [c for c in df.columns if 'PC' in c]
    for col in pc_cols:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    
    return df.dropna(subset=['PC1']), eigenval

def plot_soybean_pca(df, eigenval, output='soybean_pca_clean.png'):
    # 计算解释方差百分比
    total_ev = eigenval['val'].sum()
    ev_per = (eigenval['val'] / total_ev * 100).round(2)

    plt.figure(figsize=(10, 8))
    
    # 使用 Seaborn 绘制，cmap 选一个冷淡高级的颜色
    sns.set_style("whitegrid")
    scatter = plt.scatter(df['PC1'], df['PC2'], 
                         c=df['PC1'], cmap='coolwarm', 
                         alpha=0.7, s=60, edgecolors='w', linewidth=0.5)
    
    plt.xlabel(f'PC1 ({ev_per[0]}%)', fontsize=12)
    plt.ylabel(f'PC2 ({ev_per[1]}%)', fontsize=12)
    plt.title('Soybean Population Structure (PCA)', fontsize=15, pad=20)
    
    # 添加方差解释比例的标注
    info_text = f"PC1+PC2: {ev_per[0]+ev_per[1]:.1f}%\nN_samples: {len(df)}"
    plt.text(0.95, 0.05, info_text, transform=plt.gca().transAxes, 
             fontsize=10, verticalalignment='bottom', horizontalalignment='right',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

    plt.colorbar(scatter, label='PC1 Gradient')
    plt.savefig(output, dpi=300, bbox_inches='tight')
    print(f"✓ 图表已生成: {output}")

# 执行部分
if __name__ == "__main__":
    # 确保文件名与你服务器上的匹配
    file_prefix = 'soybean_pca' 
    eigenvec_df, eigenval_df = load_pca_data(file_prefix)
    plot_soybean_pca(eigenvec_df, eigenval_df)