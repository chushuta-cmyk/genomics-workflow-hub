import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# 读取你之前算好的均值表型文件
df = pd.read_csv('phenotype_long.tsv', sep='\t')
# 只要均值
avg = df.groupby('ID')[['Protein', 'Oil']].mean().reset_index()

# 绘制 Protein 和 Oil 的相关性散点图
plt.figure(figsize=(8, 6))
plt.scatter(avg['Protein'], avg['Oil'], alpha=0.5, c='forestgreen')
plt.xlabel('Protein Content (%)')
plt.ylabel('Oil Content (%)')
plt.title('Correlation between Protein and Oil')
plt.grid(True, linestyle='--')
plt.savefig('protein_vs_oil_scatter.png')

print("性状分布图已生成: protein_vs_oil_scatter.png")