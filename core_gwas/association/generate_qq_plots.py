import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import os

def get_p_values(file_path):
    """大文件分块读取 P 值，防止内存崩溃"""
    p_list = []
    print(f"正在读取: {os.path.basename(file_path)}")
    # 假设你的列名是 p_wald，如果不是请根据实际情况修改
    for chunk in pd.read_csv(file_path, sep='\t', usecols=['p_wald'], chunksize=200000):
        p_list.extend(chunk['p_wald'].tolist())
    return np.array(p_list)

traits = {
    'Oil': 'gwas_oil_results.clean.txt',
    'Protein': 'gwas_protein_results.clean.txt',
    'Size': 'gwas_size_results.clean.txt'
}

plt.figure(figsize=(18, 5))

for i, (name, path) in enumerate(traits.items()):
    p_vals = get_p_values(path)
    p_vals = p_vals[~np.isnan(p_vals)] # 剔除空值
    p_vals.sort()
    
    n = len(p_vals)
    expected = np.arange(1, n + 1) / (n + 1)
    
    # 转换为 -log10
    observed_log = -np.log10(p_vals)
    expected_log = -np.log10(expected[::-1]) # 降序排列以对应
    
    # 计算膨胀因子 Lambda
    # Lambda = median(chi2_observed) / median(chi2_expected)
    chi2_observed = stats.chi2.ppf(1 - p_vals, 1)
    lambda_gc = np.median(chi2_observed) / stats.chi2.ppf(0.5, 1)
    
    ax = plt.subplot(1, 3, i + 1)
    ax.scatter(expected_log, observed_log, c='royalblue', s=10, alpha=0.5, edgecolors='none')
    ax.plot([0, max(expected_log)], [0, max(expected_log)], color='red', linestyle='--')
    
    ax.set_title(f'QQ Plot: {name}\n$\lambda_{{GC}} = {lambda_gc:.3f}$')
    ax.set_xlabel('Expected $-log_{10}(P)$')
    ax.set_ylabel('Observed $-log_{10}(P)$')
    ax.grid(True, linestyle=':', alpha=0.6)

plt.tight_layout()
plt.savefig('combined_qq_plots.png', dpi=300)
print("所有 QQ 图已保存至 combined_qq_plots.png")