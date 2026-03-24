#!/usr/bin/env python3
import sys
sys.path.insert(0, '.')
from final_analysis import load_data, classify_loci, create_final_summary_table

print('Loading data...')
data = load_data()
print('Data keys:', list(data.keys()))
print('Coloc results columns:', data['coloc_results'].columns.tolist())
print('Coloc results dtypes:', data['coloc_results'].dtypes)

print('\nClassifying loci...')
classified_loci = classify_loci(data)
print('Classified loci shape:', classified_loci.shape)

print('\nCreating final summary table...')
summary_df = create_final_summary_table(data, classified_loci)
print('Summary table shape:', summary_df.shape)
print('Summary columns:', summary_df.columns.tolist())
# Check for colocalized columns
coloc_cols = [c for c in summary_df.columns if 'colocalized' in c]
print('Colocalized columns:', coloc_cols)
if coloc_cols:
    print('Sample values:', summary_df[coloc_cols].head())
    print('Dtypes:', summary_df[coloc_cols].dtypes)
print('\nDone')