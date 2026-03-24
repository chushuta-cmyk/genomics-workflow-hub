import pandas as pd
import numpy as np
# 读取清理后的结果
df = pd.read_csv('gwas_protein_results.clean.txt', sep='\t')

# 统计有效SNP数（考虑LD）
# 保守估计: 通常有效SNP数 = 实际SNP数 / 2-3（由于LD）
n_snps = len(df)
n_effective = n_snps // 2  # 保守估计

# 计算Bonferroni阈值
bonferroni_p = 0.05 / n_effective
bonferroni_log10p = -np.log10(bonferroni_p)

print(f"有效SNP数: {n_effective}")
print(f"Bonferroni P值阈值: {bonferroni_p:.2e}")
print(f"-log10(P) 阈值: {bonferroni_log10p:.2f}")