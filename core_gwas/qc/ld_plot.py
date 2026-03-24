import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# 1. 读取 LD 矩阵
# plink --r2 square 生成的文件通常是以空格分隔的
print("正在加载 LD 矩阵...")
ld_matrix = pd.read_csv('04_gemma/final_ld_matrix.ld', sep='\s+', header=None)

# 2. 绘图设置
plt.figure(figsize=(12, 10))

# 使用 Reds 色板，数值范围 0 到 1 (R² 值的范围)
# annot=False 因为 SNP 较多，标数字会乱
sns.heatmap(ld_matrix, 
            cmap='Reds', 
            vmin=0, vmax=1, 
            cbar_kws={'label': 'Linkage Disequilibrium (R²)'},
            xticklabels=False, yticklabels=False)

plt.title('Local LD Structure: Chr10 (14.5Mb - 15.5Mb)', fontsize=15, fontweight='bold')
plt.xlabel(f'SNPs within 1Mb region (n={ld_matrix.shape[0]})')
plt.ylabel('SNPs')

# 3. 保存
output_file = 'Soybean_Chr10_LD_Heatmap.png'
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"✨ 绘制完成！请查看文件: {output_file}")