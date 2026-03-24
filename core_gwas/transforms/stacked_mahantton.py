import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

traits = ['size', 'protein', 'oil']
# 使用更具质感的调色盘：深海蓝、翡翠绿、勃艮第红
colors = [['#1B2631', '#5D6D7E'], ['#0E6251', '#73C6B6'], ['#641E16', '#C0392B']] 
files = [f'output/{t}_analysis.assoc.txt' for t in traits]

fig, axes = plt.subplots(len(traits), 1, figsize=(16, 14), sharex=True, facecolor='white')

for i, (trait, file) in enumerate(zip(traits, files)):
    print(f"正在优化绘制 {trait} ...")
    df = pd.read_csv(file, sep='\t', dtype={'chr': str})
    df['chr_num'] = df['chr'].str.extract('(\d+)').astype(float)
    df = df[df['chr_num'].between(1, 20)].dropna(subset=['p_wald'])
    df['logP'] = -np.log10(df['p_wald'])
    
    # 坐标偏移计算
    df['pos_cum'] = 0
    pos_offset = 0
    ticks = []
    for chrom in range(1, 21):
        mask = (df['chr_num'] == chrom)
        if not df[mask].empty:
            df.loc[mask, 'pos_cum'] = df.loc[mask, 'ps'] + pos_offset
            ticks.append(pos_offset + (df.loc[mask, 'ps'].max() / 2))
            pos_offset += df.loc[mask, 'ps'].max()

    ax = axes[i]
    # 使用深浅交替色绘制
    for chrom in range(1, 21):
        c_df = df[df['chr_num'] == chrom]
        color = colors[i][0] if chrom % 2 != 0 else colors[i][1]
        ax.scatter(c_df['pos_cum'], c_df['logP'], c=color, s=5, alpha=0.6, edgecolors='none')
    
    # 额外叠加极显著位点 (突出显示)
    sig_df = df[df['p_wald'] < 5e-8]
    ax.scatter(sig_df['pos_cum'], sig_df['logP'], c='gold', s=12, edgecolors='black', linewidths=0.5, label='Significant')

    # 美化线条
    ax.axhline(y=-np.log10(5e-8), color='#CB4335', linestyle='-', linewidth=1.2, alpha=0.8)
    ax.axhline(y=-np.log10(1e-5), color='#F39C12', linestyle=':', linewidth=1, alpha=0.6)
    
    ax.set_ylabel(f'{trait.upper()}\n-log10(P)', fontsize=13, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', linestyle='--', alpha=0.2)

plt.xticks(ticks, [str(i) for i in range(1, 21)])
plt.xlabel('Soybean Chromosomes (1-20)', fontsize=15, fontweight='bold')
plt.suptitle('High-Resolution GWAS Landscape of Soybean Quality Traits', fontsize=20, fontweight='bold', y=0.98)

plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig('Soybean_GWAS_Visual_Masterpiece.png', dpi=400) # 调高 DPI
print("✅ 绝美图表已生成: Soybean_GWAS_Visual_Masterpiece.png")