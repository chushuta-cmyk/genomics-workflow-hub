import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
import pandas as pd

# 读取表型数据，跨年份求均值（消除年份效应）
df = pd.read_csv('phenotype_long.tsv', sep='\t')
pheno_avg = df.groupby('ID')[['Protein', 'Oil', '100SW']].mean().reset_index()
# 1. 计算相关系数矩阵 (Pearson Correlation)
corr_matrix = pheno_avg[['Protein', 'Oil', '100SW']].corr()
print("═══ 性状相关性矩阵 ═══")
print(corr_matrix)

# 2. 绘制多性状联合分布图 (Pairplot)
# 该图同时包含：对角线(分布直方图)、下方(散点图+拟合线)、上方(相关系数)
def corrfunc(x, y, **kws):
    r, p = stats.pearsonr(x, y)
    ax = plt.gca()
    ax.annotate(f"r = {r:.2f}\np = {p:.2e}",
                xy=(0.1, 0.9), xycoords=ax.transAxes,
                fontsize=10, fontweight='bold', color='red')

sns.set_style("white")
g = sns.PairGrid(pheno_avg[['Protein', 'Oil', '100SW']].dropna(), aspect=1.2)
g.map_upper(corrfunc)           # 右上方显示相关系数
g.map_lower(sns.regplot,        # 左下方显示散点图+线性回归拟合
            scatter_kws={'alpha':0.5, 's':20, 'color':'teal'},
            line_kws={'color':'red', 'lw':1.5})
g.map_diag(sns.histplot,        # 对角线显示各性状频率分布
           kde=True, color='steelblue', edgecolor='w')

plt.subplots_adjust(top=0.9)
g.fig.suptitle('Soybean Traits Relationship & Distribution Analysis', fontsize=16)

# 保存文件
plt.savefig('soybean_traits_correlation.png', dpi=300, bbox_inches='tight')
print("\n✓ 相关性综合图已生成: soybean_traits_correlation.png")

# 3. 重点检查：Protein vs Oil 散点图
plt.figure(figsize=(8, 6))
sns.scatterplot(data=pheno_avg, x='Protein', y='Oil', size='100SW', 
                hue='100SW', palette='viridis', alpha=0.7)
plt.title('Protein vs Oil: The Carbon-Nitrogen Trade-off', fontsize=14)
plt.grid(True, linestyle='--', alpha=0.6)
plt.savefig('protein_vs_oil_scatter.jpg', dpi=300)

