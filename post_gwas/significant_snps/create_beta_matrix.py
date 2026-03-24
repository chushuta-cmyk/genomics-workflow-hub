import pandas as pd
df = pd.read_csv('cross_trait_snp_effects.txt', sep='\t')
# Select columns of interest
simple = df[['rs', 'chr', 'pos', 'trait_significant', 'beta_size', 'beta_protein', 'beta_oil']]
simple.to_csv('cross_trait_beta_matrix.txt', sep='\t', index=False)
print("Beta matrix saved.")