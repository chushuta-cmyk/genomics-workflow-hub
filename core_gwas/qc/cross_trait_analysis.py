# Generalized from soybean project-specific path layout.
#!/usr/bin/env python3
import pandas as pd
import os
import sys

# Paths to GWAS files
SIZE_FILE = 'data/input/workflow/output/size_analysis.assoc.txt'
PROTEIN_FILE = 'data/input/workflow/output/protein_analysis.assoc.txt'
OIL_FILE = 'data/input/workflow/output/oil_analysis.assoc.txt'

# Output directory
OUT_DIR = '.'

def load_gwas_dict(file_path):
    """Load GWAS file into dict keyed by rs with beta and p_wald."""
    print(f"Loading {file_path}...")
    df = pd.read_csv(file_path, sep='\t')
    # Ensure rs column is string
    df['rs'] = df['rs'].astype(str)
    # Create dict: rs -> (beta, p_wald)
    gwas_dict = {}
    for _, row in df.iterrows():
        gwas_dict[row['rs']] = (row['beta'], row['p_wald'])
    print(f"Loaded {len(gwas_dict)} SNPs")
    return gwas_dict

def load_significant(file_path, trait):
    """Load significant SNPs file and add trait column."""
    df = pd.read_csv(file_path, sep='\t')
    df['trait'] = trait
    return df

def main():
    # Load significant SNPs for each trait
    size_sig = load_significant(os.path.join(OUT_DIR, 'size_significant.txt'), '100SW')
    protein_sig = load_significant(os.path.join(OUT_DIR, 'protein_significant.txt'), 'protein')
    oil_sig = load_significant(os.path.join(OUT_DIR, 'oil_significant.txt'), 'oil')
    
    # Merge all significant SNPs
    merged = pd.concat([size_sig, protein_sig, oil_sig], ignore_index=True)
    merged.to_csv(os.path.join(OUT_DIR, 'all_significant_snps.txt'), sep='\t', index=False)
    print(f"Merged {len(merged)} significant SNPs from all traits")
    
    # Load full GWAS data for cross-trait lookup
    size_dict = load_gwas_dict(SIZE_FILE)
    protein_dict = load_gwas_dict(PROTEIN_FILE)
    oil_dict = load_gwas_dict(OIL_FILE)
    
    # Prepare cross-trait table
    # For each unique SNP across all significant lists
    all_snps = set(merged['rs'].unique())
    print(f"Unique significant SNPs across traits: {len(all_snps)}")
    
    cross_data = []
    for rs in all_snps:
        # Get trait(s) where this SNP is significant
        traits = merged.loc[merged['rs'] == rs, 'trait'].tolist()
        trait_str = ','.join(traits)
        
        # Retrieve beta and p_wald from each GWAS
        size_beta, size_p = size_dict.get(rs, (None, None))
        protein_beta, protein_p = protein_dict.get(rs, (None, None))
        oil_beta, oil_p = oil_dict.get(rs, (None, None))
        
        # Also get chromosome and position from any of the files (first occurrence)
        # Use merged row
        row = merged.loc[merged['rs'] == rs].iloc[0]
        chrom = row['chr']
        pos = row['ps']
        allele1 = row['allele1']
        allele0 = row['allele0']
        af = row['af']
        
        cross_data.append({
            'rs': rs,
            'chr': chrom,
            'pos': pos,
            'allele1': allele1,
            'allele0': allele0,
            'af': af,
            'trait_significant': trait_str,
            'beta_size': size_beta,
            'p_wald_size': size_p,
            'beta_protein': protein_beta,
            'p_wald_protein': protein_p,
            'beta_oil': oil_beta,
            'p_wald_oil': oil_p
        })
    
    cross_df = pd.DataFrame(cross_data)
    cross_df.to_csv(os.path.join(OUT_DIR, 'cross_trait_snp_effects.txt'), sep='\t', index=False)
    print(f"Cross-trait SNP effect table saved with {len(cross_df)} SNPs")
    
    # Summary statistics
    print("\n=== Summary ===")
    print(f"Size significant SNPs: {len(size_sig)}")
    print(f"Protein significant SNPs: {len(protein_sig)}")
    print(f"Oil significant SNPs: {len(oil_sig)}")
    print(f"Unique significant SNPs across traits: {len(all_snps)}")
    
    # SNPs significant in multiple traits
    multi_trait = [rs for rs in all_snps if len(merged.loc[merged['rs'] == rs, 'trait'].unique()) > 1]
    print(f"SNPs significant in multiple traits: {len(multi_trait)}")
    if multi_trait:
        print("Multi-trait SNPs:", multi_trait)
    
    # Save multi-trait SNPs
    multi_df = cross_df[cross_df['rs'].isin(multi_trait)]
    multi_df.to_csv(os.path.join(OUT_DIR, 'multi_trait_significant_snps.txt'), sep='\t', index=False)

if __name__ == '__main__':
    main()