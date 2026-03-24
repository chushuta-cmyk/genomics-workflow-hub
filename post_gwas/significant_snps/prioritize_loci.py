#!/usr/bin/env python3
"""
Prioritize loci for Bayesian colocalization based on cross-trait SNP effects.

Input: cross_trait_snp_effects.txt
Output: prioritized_loci.tsv
"""

import pandas as pd
import numpy as np
import sys
from pathlib import Path

def main():
    # Paths
    input_file = Path("cross_trait_snp_effects.txt")
    output_file = Path("coloc/prioritized_loci.tsv")
    
    # Load data
    df = pd.read_csv(input_file, sep='\t')
    print(f"Loaded {len(df)} SNPs")
    
    # Ensure chromosome is string
    df['chr'] = df['chr'].astype(str)
    
    # Determine most significant p-value for each SNP across traits
    p_cols = ['p_wald_size', 'p_wald_protein', 'p_wald_oil']
    df['min_p'] = df[p_cols].min(axis=1)
    
    # Sort by most significant (smallest p-value)
    df_sorted = df.sort_values('min_p').reset_index(drop=True)
    
    # Define window size
    window = 250000  # ±250 kb
    
    # Initialize list of loci
    loci = []
    
    # Process SNPs in order of significance
    for _, row in df_sorted.iterrows():
        chrom = row['chr']
        pos = row['pos']
        start = max(0, pos - window)
        end = pos + window
        
        # Check if SNP falls within any existing locus
        found = False
        for locus in loci:
            if locus['chr'] == chrom:
                # Check overlap: new start <= locus_end AND new end >= locus_start
                if start <= locus['end'] and end >= locus['start']:
                    # Merge: expand locus boundaries
                    locus['start'] = min(locus['start'], start)
                    locus['end'] = max(locus['end'], end)
                    # Add SNP to locus SNPs list
                    locus['snps'].append(row['rs'])
                    # Update lead SNP if this SNP has lower p-value
                    if row['min_p'] < locus['lead_p']:
                        locus['lead_snp'] = row['rs']
                        locus['lead_pos'] = pos
                        locus['lead_p'] = row['min_p']
                    found = True
                    break
        
        # If not within any existing locus, create new locus
        if not found:
            # Determine which traits are significant (p < 5e-8) for this SNP
            sig_traits = []
            if row['p_wald_size'] < 5e-8:
                sig_traits.append('size')
            if row['p_wald_protein'] < 5e-8:
                sig_traits.append('protein')
            if row['p_wald_oil'] < 5e-8:
                sig_traits.append('oil')
            
            loci.append({
                'chr': chrom,
                'start': start,
                'end': end,
                'lead_snp': row['rs'],
                'lead_pos': pos,
                'lead_p': row['min_p'],
                'sig_traits': ','.join(sig_traits) if sig_traits else 'none',
                'snps': [row['rs']]
            })
    
    # Create output DataFrame
    loci_df = pd.DataFrame([
        {
            'locus_id': f"locus_{i+1:03d}",
            'chr': locus['chr'],
            'start': locus['start'],
            'end': locus['end'],
            'locus_size': locus['end'] - locus['start'],
            'lead_snp': locus['lead_snp'],
            'lead_pos': locus['lead_pos'],
            'lead_p': locus['lead_p'],
            'sig_traits': locus['sig_traits'],
            'num_snps': len(locus['snps']),
            'snps': ';'.join(locus['snps'][:10])  # First 10 SNPs only
        }
        for i, locus in enumerate(loci)
    ])
    
    # Sort loci by chromosome and start position
    loci_df = loci_df.sort_values(['chr', 'start']).reset_index(drop=True)
    
    # Save to file
    output_file.parent.mkdir(exist_ok=True)
    loci_df.to_csv(output_file, sep='\t', index=False)
    print(f"Saved {len(loci_df)} prioritized loci to {output_file}")
    
    # Print summary
    print("\nPrioritized Loci Summary:")
    print(f"  Total loci: {len(loci_df)}")
    print(f"  Loci with multiple SNPs: {(loci_df['num_snps'] > 1).sum()}")
    print(f"  Total SNPs across loci: {loci_df['num_snps'].sum()}")
    
    # Count loci by chromosome
    chrom_counts = loci_df['chr'].value_counts()
    print("\nLoci by chromosome:")
    for chrom, count in chrom_counts.items():
        print(f"  {chrom}: {count} loci")

if __name__ == "__main__":
    main()