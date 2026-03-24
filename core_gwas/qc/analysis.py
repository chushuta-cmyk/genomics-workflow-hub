# Generalized from soybean project-specific path layout.
#!/usr/bin/env python3
"""
Proteomic Validation Analysis Script
Goal: Validate EBV-driven biphasic immune response using proteomics
Author: Senior Bioinformatics Scientist
Date: 2025
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from matplotlib import rcParams
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# Set style for SCI publications
sns.set_style('whitegrid')
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
rcParams['font.size'] = 10
rcParams['axes.labelsize'] = 10
rcParams['axes.titlesize'] = 12
rcParams['xtick.labelsize'] = 9
rcParams['ytick.labelsize'] = 9
rcParams['legend.fontsize'] = 9
rcParams['figure.dpi'] = 300
rcParams['savefig.dpi'] = 300
rcParams['savefig.bbox'] = 'tight'

# Paths
workspace = 'docs/proteomic'
output_dir = os.path.join(workspace, 'pics')
os.makedirs(output_dir, exist_ok=True)

diff_file = os.path.join(workspace, '2_Diff_ResultSummary.xlsx')
quant_file = os.path.join(workspace, 'IdentficationQuantificantion_ResultSummary.xlsx')
dia_file = os.path.join(workspace, 'DIA.txt')

# Load data
print("Loading data...")
diff_df = pd.read_excel(diff_file, sheet_name=0)
quant_df = pd.read_excel(quant_file, sheet_name=0)
dia_genes = pd.read_csv(dia_file, sep='\t', header=None).iloc[:,0].tolist()

# Data Pre-processing: Map columns for Log2FC and P-values
print("\n=== Data Pre-processing ===")
print("Columns in differential result file:")
print(diff_df.columns.tolist())

# Determine Log2FC column
# We have ratio column 'NPC43-Re/NPC43-Ctl' (fold change)
fc_col = 'NPC43-Re/NPC43-Ctl'
if fc_col in diff_df.columns:
    diff_df['log2FC'] = np.log2(diff_df[fc_col])
    print(f"Computed log2FC from column '{fc_col}'")
else:
    # Look for column containing 'log2' or 'FC'
    fc_candidates = [c for c in diff_df.columns if 'log2' in c.lower() or 'fc' in c.lower()]
    if fc_candidates:
        fc_col = fc_candidates[0]
        diff_df['log2FC'] = diff_df[fc_col]
        print(f"Using column '{fc_col}' as log2FC")
    else:
        raise ValueError("Cannot find fold change column")

# Determine P-value column
pval_col = None
for col in ['Qvalue', 'Pvalue', 'pvalue', 'adj.P.Val', 'FDR']:
    if col in diff_df.columns:
        pval_col = col
        break
if pval_col is None:
    # Look for column containing 'p' or 'q'
    p_candidates = [c for c in diff_df.columns if 'p' in c.lower() or 'q' in c.lower()]
    if p_candidates:
        pval_col = p_candidates[0]
    else:
        raise ValueError("Cannot find p-value column")
print(f"Using column '{pval_col}' as p-value")

# Add -log10(p)
diff_df['neg_log10_p'] = -np.log10(diff_df[pval_col])

# Report mapping results
print(f"\nMapping results:")
print(f"- Log2FC derived from: {fc_col}")
print(f"- P-value column used: {pval_col}")
print(f"- Number of proteins: {len(diff_df)}")
print(f"- Range of log2FC: [{diff_df['log2FC'].min():.2f}, {diff_df['log2FC'].max():.2f}]")
print(f"- Range of -log10(p): [{diff_df['neg_log10_p'].min():.2f}, {diff_df['neg_log10_p'].max():.2f}]")
print(f"- Significant proteins (p < 0.05): {(diff_df[pval_col] < 0.05).sum()}")
print(f"- Proteins with |log2FC| >= 1: {(abs(diff_df['log2FC']) >= 1).sum()}")

# Save mapping report
with open(os.path.join(output_dir, 'mapping_report.txt'), 'w') as f:
    f.write("Proteomics Data Mapping Report\n")
    f.write("==============================\n")
    f.write(f"Log2FC source column: {fc_col}\n")
    f.write(f"P-value source column: {pval_col}\n")
    f.write(f"Total proteins: {len(diff_df)}\n")
    f.write(f"Significant (p<0.05): {(diff_df[pval_col] < 0.05).sum()}\n")
    f.write(f"|log2FC|>=1: {(abs(diff_df['log2FC']) >= 1).sum()}\n")
    f.write(f"Up-regulated (log2FC >=1): {(diff_df['log2FC'] >= 1).sum()}\n")
    f.write(f"Down-regulated (log2FC <= -1): {(diff_df['log2FC'] <= -1).sum()}\n")

# Volcano Plotting
print("\n=== Generating Volcano Plot ===")
# Thresholds
fc_thresh = 1
p_thresh = 0.05
neg_log10_p_thresh = -np.log10(p_thresh)

# Determine up/down significant
diff_df['significant'] = 'ns'
diff_df.loc[(diff_df[pval_col] < p_thresh) & (diff_df['log2FC'] >= fc_thresh), 'significant'] = 'up'
diff_df.loc[(diff_df[pval_col] < p_thresh) & (diff_df['log2FC'] <= -fc_thresh), 'significant'] = 'down'

# Colors
palette = {'up': '#D62728', 'down': '#1F77B4', 'ns': '#7F7F7F'}

# Genes to annotate
genes_to_annotate = ['THBS1', 'EPHA2', 'STAT3', 'S100A10', 'ITGA3']
# Find indices
annot_indices = {}
for gene in genes_to_annotate:
    matches = diff_df[diff_df['Gene'] == gene]
    if matches.empty:
        matches = diff_df[diff_df['Gene'].str.lower() == gene.lower()]
    if not matches.empty:
        idx = matches.index[0]
        annot_indices[gene] = idx
        print(f"Found {gene}: log2FC={diff_df.loc[idx, 'log2FC']:.2f}, p={diff_df.loc[idx, pval_col]:.2e}")

# Create figure for volcano only (for panel A)
fig_volc, ax_volc = plt.subplots(figsize=(6, 5))
# Scatter
for sig in ['ns', 'up', 'down']:
    subset = diff_df[diff_df['significant'] == sig]
    ax_volc.scatter(subset['log2FC'], subset['neg_log10_p'], 
                    c=palette[sig], label=sig, s=20, alpha=0.7, edgecolors='none')

# Threshold lines
ax_volc.axvline(x=fc_thresh, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
ax_volc.axvline(x=-fc_thresh, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
ax_volc.axhline(y=neg_log10_p_thresh, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)

# Annotate genes
for gene, idx in annot_indices.items():
    ax_volc.annotate(gene, 
                     xy=(diff_df.loc[idx, 'log2FC'], diff_df.loc[idx, 'neg_log10_p']),
                     xytext=(5, 5), textcoords='offset points',
                     fontsize=9, color='black',
                     arrowprops=dict(arrowstyle='-', color='black', alpha=0.5, lw=0.5))

# Labels
ax_volc.set_xlabel('Log$_2$ Fold Change (High_CP / Low_CP)', fontsize=11)
ax_volc.set_ylabel('-Log$_{10}$ Q-value', fontsize=11)
ax_volc.set_title('Differential Protein Expression', fontsize=12, pad=10)
ax_volc.legend(title='Significance', loc='upper right')
ax_volc.grid(True, linestyle=':', alpha=0.5)

# Save volcano alone
volc_path = os.path.join(output_dir, 'volcano_plot.png')
fig_volc.savefig(volc_path, dpi=300)
print(f"Saved volcano plot to {volc_path}")

# Heatmap Generation
print("\n=== Generating Heatmap ===")
# Use differentially expressed proteins (significant up/down)
sig_proteins = diff_df[diff_df['significant'].isin(['up', 'down'])]
print(f"Number of differentially expressed proteins: {len(sig_proteins)}")

# Extract abundance values (raw) from diff_df (columns NPC43-Ctl, NPC43-Re)
abundance_cols = ['NPC43-Ctl', 'NPC43-Re']
if all(col in diff_df.columns for col in abundance_cols):
    abundance = sig_proteins[abundance_cols].copy()
else:
    # Fallback to quant_df: merge by Accession or Gene
    abundance = pd.merge(sig_proteins[['Gene', 'Accession']], 
                         quant_df[['Gene', 'NPC43-Ctl', 'NPC43-Re']], 
                         on='Gene', how='left')
    abundance_cols = ['NPC43-Ctl', 'NPC43-Re']

# Ensure numeric
abundance[abundance_cols] = abundance[abundance_cols].apply(pd.to_numeric, errors='coerce')
# Replace zeros with small value to avoid log issues
abundance[abundance_cols] = abundance[abundance_cols].replace(0, np.nan)
# Log2 transform abundance? Usually heatmap uses Z-score of raw abundance; we'll Z-score log2 transformed.
# Log2 transform after adding pseudocount
log_abundance = np.log2(abundance[abundance_cols] + 1)

# Z-score normalization per protein (row-wise)
from scipy.stats import zscore
z_matrix = log_abundance.apply(lambda row: zscore(row, nan_policy='omit'), axis=1, result_type='expand')
z_matrix.columns = abundance_cols
z_matrix.index = sig_proteins['Gene'].values

# Order rows by log2FC descending (up-regulated first)
sig_proteins_sorted = sig_proteins.sort_values('log2FC', ascending=False)
z_matrix = z_matrix.loc[sig_proteins_sorted['Gene']]

# Prepare column labels (High_CP vs Low_CP)
column_labels = ['Low_CP', 'High_CP']  # NPC43-Ctl = Low, NPC43-Re = High
z_matrix.columns = column_labels

# Create heatmap figure
fig_heat, ax_heat = plt.subplots(figsize=(4, 8))
# Use seaborn heatmap
sns.heatmap(z_matrix, ax=ax_heat, cmap='RdBu_r', center=0,
            cbar_kws={'label': 'Z-score', 'shrink': 0.8},
            yticklabels=False)  # hide protein names for clarity
ax_heat.set_xlabel('EBV Load', fontsize=11)
ax_heat.set_ylabel('Differentially Expressed Proteins', fontsize=11)
ax_heat.set_title('Protein Abundance Z-scores', fontsize=12, pad=10)
# Adjust layout
fig_heat.tight_layout()
heat_path = os.path.join(output_dir, 'heatmap.png')
fig_heat.savefig(heat_path, dpi=300)
print(f"Saved heatmap to {heat_path}")

# SCI Multipanel Integration
print("\n=== Generating Multipanel Figure ===")
fig, axes = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={'width_ratios': [1, 0.6]})

# Panel A: Volcano
axA = axes[0]
for sig in ['ns', 'up', 'down']:
    subset = diff_df[diff_df['significant'] == sig]
    axA.scatter(subset['log2FC'], subset['neg_log10_p'], 
                c=palette[sig], label=sig, s=20, alpha=0.7, edgecolors='none')
axA.axvline(x=fc_thresh, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
axA.axvline(x=-fc_thresh, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
axA.axhline(y=neg_log10_p_thresh, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
for gene, idx in annot_indices.items():
    axA.annotate(gene, 
                 xy=(diff_df.loc[idx, 'log2FC'], diff_df.loc[idx, 'neg_log10_p']),
                 xytext=(5, 5), textcoords='offset points',
                 fontsize=9, color='black',
                 arrowprops=dict(arrowstyle='-', color='black', alpha=0.5, lw=0.5))
axA.set_xlabel('Log$_2$ Fold Change (High_CP / Low_CP)', fontsize=11)
axA.set_ylabel('-Log$_{10}$ Q-value', fontsize=11)
axA.set_title('A. Differential Protein Expression', fontsize=12, loc='left')
axA.legend(title='Significance', loc='upper right')
axA.grid(True, linestyle=':', alpha=0.5)

# Panel B: Heatmap
axB = axes[1]
# Plot heatmap using imshow for better control
im = axB.imshow(z_matrix, aspect='auto', cmap='RdBu_r', interpolation='nearest')
axB.set_xticks(range(len(column_labels)))
axB.set_xticklabels(column_labels, fontsize=10)
axB.set_yticks([])
axB.set_ylabel('')
axB.set_title('B. Protein Abundance Z-scores', fontsize=12, loc='left')
# Add colorbar
cbar = fig.colorbar(im, ax=axB, shrink=0.8)
cbar.set_label('Z-score', fontsize=10)

# Adjust layout
fig.tight_layout()
multipanel_path = os.path.join(output_dir, 'multipanel_figure.png')
fig.savefig(multipanel_path, dpi=300)
print(f"Saved multipanel figure to {multipanel_path}")

# Generate Figure Legend (caption)
legend_text = """
Figure 1. Proteomic validation of EBV-driven biphasic immune response in NPC tumors.
(A) Volcano plot displaying differential protein expression between High_CP (high EBV load) and Low_CP (low EBV load) tumors. Proteins with |log2FC| ≥ 1 and Q-value < 0.05 are considered significant (red: up-regulated, blue: down-regulated). Key immune‑related proteins (THBS1, EPHA2, STAT3, S100A10, ITGA3) are annotated. Dashed lines indicate significance thresholds.
(B) Heatmap of Z‑score normalized protein abundance for differentially expressed proteins (rows) across High_CP and Low_CP samples (columns). Proteins are ordered by decreasing log2 fold change. Color scale represents row‑wise Z‑scores, highlighting relative abundance patterns associated with EBV load.
"""
legend_path = os.path.join(output_dir, 'figure_legend.txt')
with open(legend_path, 'w') as f:
    f.write(legend_text.strip())
print(f"Saved figure legend to {legend_path}")

# Summary statistics for manuscript
print("\n=== Summary Statistics ===")
up = diff_df[diff_df['significant'] == 'up']
down = diff_df[diff_df['significant'] == 'down']
print(f"Up-regulated proteins: {len(up)}")
print(f"Down-regulated proteins: {len(down)}")
print(f"Top 5 up-regulated (by log2FC):")
print(up.nlargest(5, 'log2FC')[['Gene', 'log2FC', pval_col]])
print(f"Top 5 down-regulated (by log2FC):")
print(down.nsmallest(5, 'log2FC')[['Gene', 'log2FC', pval_col]])

# Save results table
results_table = diff_df[['Gene', 'Accession', 'log2FC', pval_col, 'neg_log10_p', 'significant']]
results_table.to_csv(os.path.join(output_dir, 'differential_expression_results.csv'), index=False)
print(f"Saved full results table to {os.path.join(output_dir, 'differential_expression_results.csv')}")

print("\n=== Analysis Complete ===")