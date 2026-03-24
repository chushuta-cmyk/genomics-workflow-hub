#!/usr/bin/env python3
"""
Merge nearby SNPs into loci using genomic window.

This script groups SNPs within a specified genomic window (±window bp)
into loci for downstream analysis. Each locus is defined by a lead SNP
(lowest p-value) and includes summary statistics for all SNPs in the locus.

Usage:
    python merge_loci.py --snp-table cross_trait_effect.tsv --window 250000 --output locus_summary.tsv

Input:
    Cross-trait SNP effect table from cross_trait_effects.py

Output:
    locus_summary.tsv - Summary of merged loci with lead SNP information
"""

import pandas as pd
import argparse
import sys
from pathlib import Path

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Merge nearby SNPs into loci using genomic window"
    )
    parser.add_argument(
        "--snp-table",
        required=True,
        help="Cross-trait SNP effect table (TSV format)"
    )
    parser.add_argument(
        "--window",
        type=int,
        default=250000,
        help="Genomic window size in base pairs (default: 250000)"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output file path for locus summary"
    )
    parser.add_argument(
        "--p-column",
        default="p_wald_size",
        help="Column to use for lead SNP selection (default: p_wald_size)"
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print detailed progress information"
    )
    
    return parser.parse_args()

def load_snp_table(filepath, verbose=False):
    """
    Load SNP table and prepare for locus merging.
    
    Args:
        filepath: Path to SNP effect table
        verbose: Print progress if True
    
    Returns:
        DataFrame with SNPs sorted by chromosome and position
    """
    if verbose:
        print(f"Loading SNP table from {filepath}...")
    
    df = pd.read_csv(filepath, sep='\t')
    
    # Required base columns
    required_base_cols = ['rs', 'chr']
    missing_base = [col for col in required_base_cols if col not in df.columns]
    if missing_base:
        raise ValueError(f"Missing required columns in {filepath}: {missing_base}")
    
    # Accept either 'pos' or 'ps' as genomic position
    if 'pos' not in df.columns:
        if 'ps' in df.columns:
            df = df.copy()
            df['pos'] = df['ps']
            if verbose:
                print("  Column 'pos' not found; using 'ps' as position column.")
        else:
            raise ValueError(f"Missing required position column in {filepath}: need 'pos' or 'ps'")
    
    # Ensure chr is string and pos is numeric
    df['chr'] = df['chr'].astype(str)
    df['pos'] = pd.to_numeric(df['pos'], errors='coerce')
    
    # Drop rows with missing positions
    n_before = len(df)
    df = df.dropna(subset=['pos'])
    n_after = len(df)
    
    if n_before != n_after and verbose:
        print(f"  Dropped {n_before - n_after} rows with missing/invalid positions")
    
    # Convert to integer if safe
    df['pos'] = df['pos'].astype(int)
    
    # Sort by chromosome and position
    df = df.sort_values(['chr', 'pos']).reset_index(drop=True)
    
    if verbose:
        print(f"  Loaded {len(df):,} SNPs")
        print(f"  Chromosomes: {df['chr'].unique()[:10]}")
        if len(df['chr'].unique()) > 10:
            print(f"  ... and {len(df['chr'].unique()) - 10} more")
    
    return df

def merge_snps_into_loci(df, window_size, p_column, verbose=False):
    """
    Merge SNPs into loci based on genomic window.
    
    Args:
        df: DataFrame with SNPs sorted by chromosome and position
        window_size: Window size in base pairs
        p_column: Column to use for lead SNP selection
        verbose: Print progress if True
    
    Returns:
        DataFrame with loci summary
    """
    if verbose:
        print(f"Merging SNPs into loci with window ±{window_size:,} bp...")
    
    loci = []
    current_locus = None
    
    for i, row in df.iterrows():
        chr_val = row['chr']
        pos_val = row['pos']
        
        # If first SNP or different chromosome, start new locus
        if current_locus is None or chr_val != current_locus['chr']:
            if current_locus is not None:
                loci.append(current_locus)
            
            current_locus = {
                'chr': chr_val,
                'start': pos_val,
                'end': pos_val,
                'snps': [row.to_dict()],
                'num_snps': 1
            }
        
        # If within window of current locus, add to it
        elif pos_val - current_locus['end'] <= window_size:
            current_locus['snps'].append(row.to_dict())
            current_locus['end'] = pos_val
            current_locus['num_snps'] += 1
        
        # Otherwise, start new locus
        else:
            loci.append(current_locus)
            current_locus = {
                'chr': chr_val,
                'start': pos_val,
                'end': pos_val,
                'snps': [row.to_dict()],
                'num_snps': 1
            }
        
        # Progress update
        if verbose and i > 0 and i % 10000 == 0:
            print(f"  Processed {i:,} SNPs, found {len(loci):,} loci...")
    
    # Don't forget the last locus
    if current_locus is not None:
        loci.append(current_locus)
    
    if verbose:
        print(f"  Total loci identified: {len(loci):,}")
    
    return loci

def summarize_loci(loci, p_column, verbose=False):
    """
    Create summary DataFrame from loci.
    
    Args:
        loci: List of locus dictionaries
        p_column: Column to use for lead SNP selection
        verbose: Print progress if True
    
    Returns:
        DataFrame with locus summary
    """
    if verbose:
        print(f"Summarizing loci...")
    
    summary_rows = []
    
    for i, locus in enumerate(loci):
        # Extract SNPs in this locus
        snps = locus['snps']
        
        # Find lead SNP (lowest p-value)
        lead_snp = None
        lead_p = float('inf')
        
        for snp in snps:
            # Try to get p-value from specified column
            if p_column in snp and pd.notna(snp[p_column]):
                p_val = snp[p_column]
                if p_val < lead_p:
                    lead_p = p_val
                    lead_snp = snp
        
        # If no p-value in specified column, use first SNP
        if lead_snp is None:
            lead_snp = snps[0]
            lead_p = None
        
        # Get traits represented in this locus
        traits = set()
        for snp in snps:
            if 'trait_significant' in snp and pd.notna(snp['trait_significant']):
                trait_list = str(snp['trait_significant']).split(',')
                traits.update([t.strip() for t in trait_list if t.strip()])
            elif 'trait' in snp and pd.notna(snp['trait']):
                traits.add(str(snp['trait']).strip())
        
        traits_str = ','.join(sorted(traits)) if traits else "none"
        
        # Get all SNP IDs in locus
        snp_ids = [snp['rs'] for snp in snps if 'rs' in snp and pd.notna(snp['rs'])]
        
        # Create summary row
        summary_rows.append({
            'locus_id': f"locus_{i+1:04d}",
            'chr': locus['chr'],
            'start': locus['start'],
            'end': locus['end'],
            'locus_size': locus['end'] - locus['start'],
            'num_snps': locus['num_snps'],
            'lead_snp': lead_snp.get('rs', 'unknown'),
            'lead_snp_pos': lead_snp.get('pos', 'unknown'),
            'lead_p': lead_p,
            'traits_in_locus': traits_str,
            'snp_list': ';'.join(snp_ids[:10]),  # First 10 SNPs only
            'total_snps_in_list': len(snp_ids)
        })
    
    # Create DataFrame
    summary_df = pd.DataFrame(summary_rows)
    
    # Add locus width category
    def categorize_width(size):
        if size == 0:
            return "single_snp"
        elif size <= 10000:
            return "narrow"
        elif size <= 100000:
            return "medium"
        else:
            return "broad"
    
    summary_df['locus_width'] = summary_df['locus_size'].apply(categorize_width)
    
    if verbose:
        print(f"  Created summary for {len(summary_df):,} loci")
    
    return summary_df

def main():
    """Main execution function."""
    args = parse_arguments()
    
    if args.verbose:
        print(f"Locus Merging Analysis")
        print(f"=====================")
        print(f"SNP table: {args.snp_table}")
        print(f"Window size: ±{args.window:,} bp")
        print(f"P-value column for lead SNP: {args.p_column}")
        print(f"Output file: {args.output}")
        print()
    
    # Load SNP table
    try:
        df = load_snp_table(args.snp_table, verbose=args.verbose)
    except Exception as e:
        print(f"Error loading SNP table: {str(e)}")
        sys.exit(1)
    
    # Check if p-column exists
    if args.p_column not in df.columns:
        print(f"Warning: P-value column '{args.p_column}' not found in SNP table")
        print(f"Available columns: {', '.join(df.columns.tolist())}")
        # Use first available p-value column
        p_cols = [col for col in df.columns if col.startswith('p_wald_')]
        if p_cols:
            args.p_column = p_cols[0]
            if args.verbose:
                print(f"Using alternative p-column: {args.p_column}")
        else:
            print("Error: No p-value columns found in SNP table")
            sys.exit(1)
    
    # Merge SNPs into loci
    try:
        loci = merge_snps_into_loci(df, args.window, args.p_column, verbose=args.verbose)
    except Exception as e:
        print(f"Error merging SNPs into loci: {str(e)}")
        sys.exit(1)
    
    # Create locus summary
    try:
        summary_df = summarize_loci(loci, args.p_column, verbose=args.verbose)
    except Exception as e:
        print(f"Error creating locus summary: {str(e)}")
        sys.exit(1)
    
    # Save output
    try:
        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        summary_df.to_csv(output_path, sep='\t', index=False)
        
        if args.verbose:
            print(f"\nSaved locus summary to: {output_path}")
    except Exception as e:
        print(f"Error saving output file: {str(e)}")
        sys.exit(1)
    
    # Generate statistics
    # Calculate statistics needed for both verbose output and stats file
    total_snps = df['rs'].nunique()
    total_loci = len(summary_df)
    snps_in_loci = summary_df['num_snps'].sum()
    multi_trait = summary_df[summary_df['traits_in_locus'].str.contains(',')]
    
    if args.verbose:
        print(f"\nLocus Summary Statistics")
        print("-" * 60)
        
        print(f"Total unique SNPs: {total_snps:,}")
        print(f"Total loci identified: {total_loci:,}")
        print(f"SNPs assigned to loci: {snps_in_loci:,} ({snps_in_loci/total_snps:.1%})")
        
        # Locus size distribution
        print(f"\nLocus size distribution:")
        size_stats = summary_df['locus_size'].describe()
        print(f"  Min: {size_stats['min']:,} bp")
        print(f"  Mean: {size_stats['mean']:,.0f} bp")
        print(f"  Median: {summary_df['locus_size'].median():,.0f} bp")
        print(f"  Max: {size_stats['max']:,} bp")
        
        # SNP count per locus
        print(f"\nSNPs per locus:")
        snp_stats = summary_df['num_snps'].describe()
        print(f"  Min: {snp_stats['min']:,} SNPs")
        print(f"  Mean: {snp_stats['mean']:.1f} SNPs")
        print(f"  Median: {summary_df['num_snps'].median():.1f} SNPs")
        print(f"  Max: {snp_stats['max']:,} SNPs")
        
        # Locus width categories
        print(f"\nLocus width categories:")
        width_counts = summary_df['locus_width'].value_counts().sort_index()
        for width, count in width_counts.items():
            percentage = 100 * count / total_loci
            print(f"  {width:12s}: {count:4d} loci ({percentage:.1f}%)")
        
        # Traits per locus
        print(f"\nTraits per locus:")
        trait_counts = {}
        for traits in summary_df['traits_in_locus']:
            if traits != 'none':
                trait_list = traits.split(',')
                for trait in trait_list:
                    trait_counts[trait] = trait_counts.get(trait, 0) + 1
        
        for trait, count in sorted(trait_counts.items()):
            percentage = 100 * count / total_loci
            print(f"  {trait:12s}: {count:4d} loci ({percentage:.1f}%)")
        
        print(f"\nMulti-trait loci: {len(multi_trait):,} ({100*len(multi_trait)/total_loci:.1f}%)")
        
        print("-" * 60)
    
    # Save detailed statistics
    stats_path = Path(args.output).with_suffix('.stats.tsv')
    
    stats_data = {
        'metric': [
            'total_snps', 'total_loci', 'snps_in_loci', 'snps_not_in_loci',
            'mean_locus_size', 'median_locus_size', 'max_locus_size',
            'mean_snps_per_locus', 'median_snps_per_locus', 'max_snps_per_locus',
            'single_snp_loci', 'narrow_loci', 'medium_loci', 'broad_loci',
            'single_trait_loci', 'multi_trait_loci'
        ],
        'value': [
            total_snps, total_loci, snps_in_loci, total_snps - snps_in_loci,
            summary_df['locus_size'].mean(), summary_df['locus_size'].median(), summary_df['locus_size'].max(),
            summary_df['num_snps'].mean(), summary_df['num_snps'].median(), summary_df['num_snps'].max(),
            len(summary_df[summary_df['locus_width'] == 'single_snp']),
            len(summary_df[summary_df['locus_width'] == 'narrow']),
            len(summary_df[summary_df['locus_width'] == 'medium']),
            len(summary_df[summary_df['locus_width'] == 'broad']),
            len(summary_df[~summary_df['traits_in_locus'].str.contains(',')]),
            len(multi_trait)
        ]
    }
    
    stats_df = pd.DataFrame(stats_data)
    stats_df.to_csv(stats_path, sep='\t', index=False)
    
    if args.verbose:
        print(f"\nDetailed statistics saved to: {stats_path}")
        print(f"\nAnalysis completed successfully")

if __name__ == "__main__":
    main()