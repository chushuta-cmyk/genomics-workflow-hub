import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

traits = ['oil', 'protein', 'size']
colors = ['#2E4053', '#1D8348', '#943126'] 
files = ['oil_smr_ready.txt', 'protein_smr_ready.txt', 'size_smr_ready.txt']

fig, axes = plt.subplots(len(traits), 1, figsize=(16, 10), sharex=True)

for i, (trait, file) in enumerate(zip(traits, files)):
    print(f"正在绘制 {trait} 的叠加图层...")
    df = pd.read_csv(file, sep='\t')
    
    # 强制转换染色体为整型，并过滤非主染色体
    df['CHR'] = pd.to_numeric(df['CHR'], errors='coerce')
    df = df[df['CHR'].between(1, 20)].dropna(subset=['CHR'])
    df['CHR'] = df['CHR'].astype(int)
    
    # 计算 -log10P
    df['logP'] = -np.log10(df['p'])
    
    # 核心修复：计算全基因组连续坐标
    df['pos_cum'] = 0
    pos_offset = 0
    ticks = []
    tick_labels = []
    
    for chrom in range(1, 21):
        mask = (df['CHR'] == chrom)
        if not df[mask].empty:
            df.loc[mask, 'pos_cum'] = df.loc[mask, 'BP'] + pos_offset
            # 记录刻度位置
            ticks.append(pos_offset + (df.loc[mask, 'BP'].max() / 2))
            tick_labels.append(str(chrom))
            pos_offset += df.loc[mask, 'BP'].max()

    ax = axes[i]
    # 隔色绘图
    for chrom in range(1, 21):
        c_df = df[df['CHR'] == chrom]
        color = colors[i] if chrom % 2 != 0 else '#BDC3C7' # 奇数用主题色，偶数用灰色
        ax.scatter(c_df['pos_cum'], c_df['logP'], c=color, s=2, alpha=0.7)
    
    ax.axhline(y=6.0, color='red', linestyle='--', linewidth=0.8) # 显著线
    ax.set_ylabel(f'{trait.upper()}\n-log10(P)')
    ax.set_ylim(0, df['logP'].max() + 1) # 动态调整高度防止压扁

# 设置横轴刻度
plt.xticks(ticks, tick_labels)
plt.xlabel('Soybean Chromosome (1-20)')
plt.tight_layout()
plt.savefig('Success_Stacked_Manhattan.png', dpi=300)
print("绘制成功！请检查 Success_Stacked_Manhattan.png")