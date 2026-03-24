# Generalized from soybean project-specific path layout.
#!/usr/bin/env python3
"""
Create fine mapping panel figure.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

base_dir = "results/cultivated/fine_mapping/log_ratio_17_37608516"
selected_chr = "17"
selected_pos = 37608516
window = 250000

# Load data
region_df = pd.read_csv(os.path.join(base_dir, "local_gwas_region.tsv"), sep='\t')
ld_matrix = pd.read_csv(os.path.join(base_dir, "local_ld_matrix.tsv"), sep='\t', header=None)
genes_df = pd.read_csv(os.path.join(base_dir, "candidate_genes_in_window.tsv"), sep='\t')

# Create figure
fig = plt.figure(figsize=(12, 10))
# Top panel: regional Manhattan
ax1 = plt.subplot2grid((10, 1), (0, 0), rowspan=5)
region_df['log10p'] = -np.log10(region_df['p_wald'])
lead_mask = region_df['ps'] == selected_pos
ax1.scatter(region_df['ps'], region_df['log10p'], s=20, alpha=0.6, c='steelblue', edgecolors='none')
if lead_mask.any():
    ax1.scatter(region_df.loc[lead_mask, 'ps'], region_df.loc[lead_mask, 'log10p'], 
                s=100, c='red', edgecolors='black', label='Lead SNP')
ax1.axhline(y=-np.log10(5e-8), color='r', linestyle='--', alpha=0.5, label='p=5e-8')
ax1.set_ylabel('-log10(p)')
ax1.set_title(f'Fine mapping of log_ratio locus chr{selected_chr}:{selected_pos}')
ax1.legend()
# Remove x-axis labels for top panel
ax1.set_xticklabels([])

# Middle panel: LD heatmap
ax2 = plt.subplot2grid((10, 1), (5, 0), rowspan=4)
sns.heatmap(ld_matrix, cmap='RdBu_r', center=0, square=True, 
            xticklabels=False, yticklabels=False, ax=ax2, cbar_kws={'label': 'r²'})
ax2.set_ylabel('SNP index')
# Add chromosome position ticks for a subset of SNPs
step = max(1, len(region_df) // 10)
tick_positions = np.arange(0, len(region_df), step)
tick_labels = region_df.iloc[tick_positions]['ps'].astype(str).tolist()
ax2.set_yticks(tick_positions)
ax2.set_yticklabels(tick_labels, fontsize=8)
ax2.set_xticks(tick_positions)
ax2.set_xticklabels(tick_labels, rotation=90, fontsize=8)
ax2.set_xlabel('Chromosome position (bp)')

# Bottom panel: gene track (simple)
ax3 = plt.subplot2grid((10, 1), (9, 0), rowspan=1)
if len(genes_df) > 0:
    for _, gene in genes_df.iterrows():
        start = gene['start']
        end = gene['end']
        strand = gene['strand']
        # draw a rectangle
        ax3.plot([start, end], [0, 0], color='black', linewidth=2)
        # arrow for strand direction
        if strand == '+':
            ax3.arrow(end, 0, -1000, 0, head_width=0.2, head_length=1000, fc='green', ec='green')
        elif strand == '-':
            ax3.arrow(start, 0, 1000, 0, head_width=0.2, head_length=1000, fc='red', ec='red')
        # gene name
        name = gene.get('gene_name', gene['gene_id'])
        ax3.text((start+end)/2, 0.3, name, ha='center', fontsize=8, rotation=90)
ax3.set_xlim(selected_pos - window, selected_pos + window)
ax3.set_ylim(-0.5, 1)
ax3.set_yticks([])
ax3.set_xlabel(f'Chromosome {selected_chr} position (bp)')
ax3.set_ylabel('Genes')

plt.tight_layout()
panel_path = os.path.join(base_dir, "fine_mapping_panel.png")
fig.savefig(panel_path, dpi=300)
print(f"Saved fine mapping panel to {panel_path}")