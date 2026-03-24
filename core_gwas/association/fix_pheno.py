#!/usr/bin/env python3
"""
Fix phenotype file format for GEMMA compatibility
GEMMA expects: FID IID phenotype1 [phenotype2 ...]
"""

import sys
import pandas as pd

def fix_phenotype_file(fam_file, pheno_file, output_file):
    """
    Ensure phenotype file matches GEMMA requirements
    """
    # Read .fam file to get correct IDs
    print("Reading .fam file...")
    fam = pd.read_csv(fam_file, sep=r'\s+', header=None,
                      names=['FID', 'IID', 'PID', 'MID', 'Sex', 'Pheno'])
    print(f"Found {len(fam)} individuals in .fam file")
    print(f"Sample IDs: {fam[['FID', 'IID']].head()}")
    
    # Read phenotype file
    print("\nReading phenotype file...")
    pheno = pd.read_csv(pheno_file, sep=r'\s+', header=None)
    print(f"Phenotype file shape: {pheno.shape}")
    print(f"First few rows:\n{pheno.head()}")
    
    # Determine phenotype file format
    if pheno.shape[1] >= 3:
        # Assume format: FID IID pheno1 [pheno2 ...]
        pheno.columns = ['FID', 'IID'] + [f'pheno{i}' for i in range(1, pheno.shape[1]-1)]
        print(f"\nAssumed columns: {list(pheno.columns)}")
    else:
        print("\nERROR: Phenotype file has fewer than 3 columns!")
        sys.exit(1)
    
    # Merge with .fam to ensure matching IDs
    print("\nMatching phenotype with .fam IDs...")
    merged = pd.merge(fam[['FID', 'IID']], pheno, on=['FID', 'IID'], how='inner')
    print(f"Matched {len(merged)} individuals")
    
    if len(merged) == 0:
        print("\nWARNING: No matching IDs found!")
        print("\nSample .fam IDs:")
        print(fam[['FID', 'IID']].head(3))
        print("\nSample phenotype IDs:")
        print(pheno[['FID', 'IID']].head(3))
        
        # Try alternative matching strategies
        print("\nTrying alternative ID matching...")
        
        # Check if IID alone matches
        fam_iid = set(fam['IID'])
        pheno_iid = set(pheno['IID'])
        iid_matches = fam_iid & pheno_iid
        print(f"IID-only matches: {len(iid_matches)}")
        
        if len(iid_matches) > 0:
            print("\nIDs match on IID but not FID. Creating compatible file...")
            # Use IID for both FID and IID in phenotype file
            pheno_fixed = pheno.copy()
            pheno_cols = [col for col in pheno.columns if col.startswith('pheno')]
            
            # Match on IID only
            merged = pd.merge(fam[['FID', 'IID']], 
                            pheno[['IID'] + pheno_cols], 
                            on='IID', how='inner')
            print(f"Matched {len(merged)} individuals using IID")
    
    if len(merged) == 0:
        print("\nERROR: Still no matches. Please check your ID formatting.")
        sys.exit(1)
    
    # Check for missing phenotypes
    pheno_cols = [col for col in merged.columns if col.startswith('pheno')]
    for col in pheno_cols:
        n_missing = merged[col].isna().sum()
        n_neg9 = (merged[col] == -9).sum()
        print(f"\n{col}: {n_missing} NA, {n_neg9} coded as -9")
    
    # Save output
    print(f"\nSaving to {output_file}...")
    merged.to_csv(output_file, sep='\t', header=False, index=False, na_rep='NA')
    print(f"Successfully wrote {len(merged)} individuals with {len(pheno_cols)} phenotype(s)")
    
    return merged

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python fix_pheno.py <fam_file> <pheno_file> <output_file>")
        print("Example: python fix_pheno.py file.fam pheno.txt output.pheno")
        sys.exit(1)
    
    fix_phenotype_file(sys.argv[1], sys.argv[2], sys.argv[3])
