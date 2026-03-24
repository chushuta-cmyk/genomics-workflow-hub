#!/usr/bin/env python3
"""
Test the bug fix for boolean coloc values.
"""
import pandas as pd
import numpy as np
import sys

# Load coloc results
coloc_df = pd.read_csv('coloc/coloc_results.tsv', sep='\t')
print('coloc_df dtypes:', coloc_df.dtypes)
print('colocalized column sample:', coloc_df['colocalized'].head())
print('colocalized dtype:', coloc_df['colocalized'].dtype)

# Simulate the problematic assignment
summary_df = pd.DataFrame({'rs': ['test1', 'test2']})
coloc_by_locus = {}
for _, row in coloc_df.iterrows():
    locus_id = row['locus_id']
    trait_pair = f"{row['trait1']}_vs_{row['trait2']}"
    if locus_id not in coloc_by_locus:
        coloc_by_locus[locus_id] = {}
    coloc_by_locus[locus_id][f'coloc_{trait_pair}_pp4'] = row['pp4']
    coloc_by_locus[locus_id][f'coloc_{trait_pair}_colocalized'] = row['colocalized']

# Add columns
coloc_columns = []
for locus_id, coloc_data in coloc_by_locus.items():
    for col in coloc_data:
        if col not in coloc_columns:
            coloc_columns.append(col)

for col in coloc_columns:
    summary_df[col] = np.nan

# Try assignment
for idx, row in summary_df.iterrows():
    locus_id = 'locus_002'  # dummy
    if locus_id in coloc_by_locus:
        for col, value in coloc_by_locus[locus_id].items():
            summary_df.at[idx, col] = value

print('summary_df dtypes:', summary_df.dtypes)
print('summary_df shape:', summary_df.shape)
print('Done')