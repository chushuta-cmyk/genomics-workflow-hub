import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

traits = ['protein', 'oil', 'size']
files = {t: 'gwas_{}_results.clean.txt'.format(t) for t in traits}

def plot_manhattan(trait, file_path):
    print("Processing {}...".format(trait))
    df = pd.read_csv(file_path, sep='\t')
    df['-log10p'] = -np.log10(df['p_wald'])
    
    # 简单的染色体坐标偏移计算，让不同染色体连在一起
    df = df.sort_values(['chr', 'ps'])
    df['ind'] = range(len(df))
    df_grouped = df.groupby('chr')['ind'].median()

    fig = plt.figure(figsize=(12, 5))
    ax = fig.add_subplot(111)
    
    colors = ['#27213C', '#94A187'] # 两种颜色交替显示染色体
    for i, (name, group) in enumerate(df.groupby('chr')):
        group.plot(kind='scatter', x='ind', y='-log10p', color=colors[i % 2], ax=ax, s=10)
    
    ax.set_xticks(df_grouped.values)
    ax.set_xticklabels(df_grouped.index)
    ax.set_xlabel('Chromosome')
    ax.set_ylabel('-log10(P)')
    ax.axhline(y=-np.log10(1e-5), color='r', linestyle='--', alpha=0.5) # 提示线
    
    plt.title('Manhattan Plot: {}'.format(trait.upper()))
    plt.savefig('{}_manhattan.png'.format(trait))
    plt.close()
    print("Saved {}_manhattan.png".format(trait))

for trait, path in files.items():
    if os.path.exists(path):
        plot_manhattan(trait, path)