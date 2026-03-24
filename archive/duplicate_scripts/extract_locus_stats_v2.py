#!/usr/bin/env python3
"""
Extract summary statistics for prioritized loci from original GWAS files.
Includes p_wald and ensures MAF is minor allele frequency.
"""

import pandas as pd
import numpy as np
import subprocess
import sys
from pathlib import Path

# Configuration
GWAS_FILES = {
    'size': '/data03/karama/soybean_gwas_filtered/workflow/output/size_analysis.assoc.txt',
    'protein': '/data03/karama/soybean_gwas_filtered/workflow/output/protein_analysis.assoc.txt',
    'oil': '/data03/karama/soybean_gwas_filtered/workflow/output/oil_analysis.assoc.txt'
}

# Sample sizes from log files
SAMPLE_SIZES = {
    'size': 374,
    'protein': 373,
    'oil': 375
}

# Column names we need
COLUMNS = ['chr', 'rs', 'ps', 'allele1', 'allele0', 'af', 'beta', 'se', 'n_miss', 'p_wald']

def extract_locus_gwas(trait, locus_chr, locus_start, locus_end):
    """Extract GWAS data for a locus using awk."""
    gwas_file = GWAS_FILES[trait]
    
    # Build awk command
    # Columns: 1=chr, 2=rs, 3=ps, 4=n_miss, 5=allele1, 6=allele0, 7=af, 8=beta, 9=se, 13=p_wald
    awk_cmd = f"""awk -F'\\t' '
BEGIN {{OFS="\\t"}}
NR==1 {{print "chr", "rs", "ps", "allele1", "allele0", "af", "beta", "se", "n_miss", "p_wald"}}
$1=="{locus_chr}" && $3>={locus_start} && $3<={locus_end} {{print $1, $2, $3, $5, $6, $7, $8, $9, $4, $13}}
' "{gwas_file}" """
    
    try:
        result = subprocess.run(awk_cmd, shell=True, capture_output=True, text=True, check=True)
        df = pd.read_csv(pd.io.common.StringIO(result.stdout), sep='\t')
        return df
    except subprocess.CalledProcessError as e:
        print(f"Error extracting locus {locus_chr}:{locus_start}-{locus_end} for {trait}: {e}")
        return pd.DataFrame(columns=COLUMNS)

def prepare_coloc_input(df, trait):
    """Prepare DataFrame for coloc input."""
    if df.empty:
        return df
    
    # Calculate sample size per SNP
    n_total = SAMPLE_SIZES[trait]
    df['N'] = n_total - df['n_miss']
    
    # Calculate variance
    df['varbeta'] = df['se'] ** 2
    
    # Ensure MAF is minor allele frequency and beta is effect of minor allele
    # af is frequency of allele1 (effect allele in GWAS)
    # If af > 0.5, allele1 is major allele, so we need to flip:
    # MAF = 1 - af, beta = -beta, p_wald unchanged, alleles swapped? We'll keep alleles as is.
    flip_mask = df['af'] > 0.5
    if flip_mask.any():
        print(f"    Flipping {flip_mask.sum()} SNPs where allele1 frequency > 0.5")
        df.loc[flip_mask, 'MAF'] = 1 - df.loc[flip_mask, 'af']
        df.loc[flip_mask, 'beta'] = -df.loc[flip_mask, 'beta']
        # Note: allele1 and allele0 remain as reported in original GWAS
        # but the effect direction is now relative to minor allele (which may be allele0)
    else:
        df['MAF'] = df['af']
    
    # Rename columns to coloc expected names
    df = df.rename(columns={
        'rs': 'snp',
        'ps': 'position',
        'p_wald': 'pvalues'
    })
    
    # Select columns needed for coloc
    coloc_cols = ['snp', 'position', 'beta', 'varbeta', 'MAF', 'N', 'pvalues']
    # Add type
    df['type'] = 'quant'
    coloc_cols.append('type')
    
    # Keep allele information for harmonization
    df['allele1'] = df['allele1']
    df['allele0'] = df['allele0']
    
    return df[coloc_cols + ['allele1', 'allele0']]

def main():
    # Paths
    loci_file = Path("coloc/prioritized_loci.tsv")
    output_dir = Path("coloc/input_v2")
    output_dir.mkdir(exist_ok=True)
    
    # Load prioritized loci
    loci_df = pd.read_csv(loci_file, sep='\t')
    print(f"Processing {len(loci_df)} loci")
    
    # Process each locus
    for _, locus in loci_df.iterrows():
        locus_id = locus['locus_id']
        chrom = str(locus['chr'])
        start = locus['start']
        end = locus['end']
        
        print(f"  Processing {locus_id}: {chrom}:{start}-{end}")
        
        # Extract data for each trait
        for trait in GWAS_FILES.keys():
            # Extract GWAS data
            gwas_df = extract_locus_gwas(trait, chrom, start, end)
            
            if gwas_df.empty:
                print(f"    Warning: No SNPs found for {trait}")
                continue
            
            # Prepare coloc input
            coloc_df = prepare_coloc_input(gwas_df, trait)
            
            # Save to file
            output_file = output_dir / f"{trait}_{locus_id}.tsv"
            coloc_df.to_csv(output_file, sep='\t', index=False)
            print(f"    Saved {len(coloc_df)} SNPs for {trait} to {output_file}")
    
    print("\nExtraction complete.")
    print(f"Output directory: {output_dir}")

if __name__ == "__main__":
    main()