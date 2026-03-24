#!/usr/bin/env python3
"""
Extract cross-trait effects for significant SNPs.

For each significant SNP (from any trait), retrieve its effect estimates
(beta, p_wald) from all GWAS results. This enables comparison of SNP effects
across multiple traits.

Usage:
    python cross_trait_effects.py --gwas-tables size.tsv protein.tsv oil.tsv --trait-names size protein oil --significant-snps sig_all_traits.tsv --output cross_trait_effect.tsv

Input:
    Simplified GWAS TSV files for all traits
    Combined significant SNPs file (sig_all_traits.tsv)

Output:
    cross_trait_effect.tsv - Table with each SNP's effects across all traits
"""

import pandas as pd
import argparse
import os
import sys
from pathlib import Path

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Extract cross-trait effects for significant SNPs"
    )
    parser.add_argument(
        "--gwas-tables",
        nargs="+",
        required=True,
        help="List of GWAS table files (TSV format) for all traits"
    )
    parser.add_argument(
        "--trait-names",
        nargs="+",
        required=True,
        help="Corresponding trait names for each GWAS table"
    )
    parser.add_argument(
        "--significant-snps",
        required=True,
        help="File with all significant SNPs (sig_all_traits.tsv)"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output file path for cross-trait effects"
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print detailed progress information"
    )
    
    return parser.parse_args()

def validate_inputs(gwas_files, trait_names, sig_snps_file):
    """Validate input files and parameters."""
    if len(gwas_files) != len(trait_names):
        raise ValueError(f"Number of GWAS files ({len(gwas_files)}) does not match number of trait names ({len(trait_names)})")
    
    for filepath in gwas_files:
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"GWAS file not found: {filepath}")
    
    if not os.path.exists(sig_snps_file):
        raise FileNotFoundError(f"Significant SNPs file not found: {sig_snps_file}")
    
    # Check for duplicate trait names
    if len(trait_names) != len(set(trait_names)):
        duplicates = [trait for trait in trait_names if trait_names.count(trait) > 1]
        raise ValueError(f"Duplicate trait names: {duplicates}")

def load_gwas_data(filepath, trait_name, verbose=False):
    """
    Load GWAS data and create lookup dictionary.
    
    Args:
        filepath: Path to GWAS TSV file
        trait_name: Name of the trait
        verbose: Print progress if True
    
    Returns:
        Dictionary {rs: (beta, p_wald, chr, pos, allele1, allele0, af)}
    """
    if verbose:
        print(f"Loading GWAS data for {trait_name}...")
    
    df = pd.read_csv(filepath, sep='\t')
    
    # Check required columns
    required_cols = ['rs', 'chr', 'ps', 'allele1', 'allele0', 'af', 'beta', 'p_wald']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns in {filepath}: {missing_cols}")
    
    # Ensure rs is string
    df['rs'] = df['rs'].astype(str)
    
    # Create dictionary
    gwas_dict = {}
    for _, row in df.iterrows():
        gwas_dict[row['rs']] = (
            row['beta'], row['p_wald'], row['chr'], row['ps'],
            row['allele1'], row['allele0'], row['af']
        )
    
    if verbose:
        print(f"  Loaded {len(gwas_dict):,} SNPs for {trait_name}")
    
    return gwas_dict

def load_significant_snps(filepath, verbose=False):
    """
    Load significant SNPs and identify which traits they're significant for.
    
    Args:
        filepath: Path to significant SNPs file
        verbose: Print progress if True
    
    Returns:
        DataFrame with significant SNPs and trait information
    """
    if verbose:
        print(f"Loading significant SNPs from {filepath}...")
    
    df = pd.read_csv(filepath, sep='\t')
    
    # Check required columns
    required_cols = ['rs', 'trait']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns in {filepath}: {missing_cols}")
    
    # Ensure rs is string
    df['rs'] = df['rs'].astype(str)
    
    if verbose:
        print(f"  Loaded {len(df):,} significant SNP-trait associations")
    
    return df

def create_cross_trait_table(sig_snps_df, gwas_dicts, trait_names, verbose=False):
    """
    Create cross-trait effect table.
    
    Args:
        sig_snps_df: DataFrame with significant SNPs
        gwas_dicts: List of GWAS dictionaries (one per trait)
        trait_names: List of trait names
        verbose: Print progress if True
    
    Returns:
        DataFrame with cross-trait effects
    """
    if verbose:
        print("Creating cross-trait effect table...")
    
    # Get unique significant SNPs across all traits
    unique_snps = sig_snps_df['rs'].unique()
    
    if verbose:
        print(f"  Unique significant SNPs: {len(unique_snps):,}")
    
    # Prepare data for output
    rows = []
    
    for rs in unique_snps:
        # Get traits where this SNP is significant
        sig_traits = sig_snps_df.loc[sig_snps_df['rs'] == rs, 'trait'].unique()
        trait_str = ','.join(sorted(sig_traits))
        
        # Initialize with None values
        data = {'rs': rs, 'trait_significant': trait_str}
        
        # Try to get chromosome, position, alleles from first GWAS file that has this SNP
        chr_val = pos_val = allele1_val = allele0_val = af_val = None
        
        for trait_name, gwas_dict in zip(trait_names, gwas_dicts):
            if rs in gwas_dict:
                beta, p_wald, chr_val, pos_val, allele1_val, allele0_val, af_val = gwas_dict[rs]
                break
        
        # If SNP not found in any GWAS dict (shouldn't happen), skip
        if chr_val is None:
            if verbose:
                print(f"  Warning: SNP {rs} not found in any GWAS data")
            continue
        
        data['chr'] = chr_val
        data['pos'] = pos_val
        data['allele1'] = allele1_val
        data['allele0'] = allele0_val
        data['af'] = af_val
        
        # Get effects from each trait
        for trait_name, gwas_dict in zip(trait_names, gwas_dicts):
            if rs in gwas_dict:
                beta, p_wald, _, _, _, _, _ = gwas_dict[rs]
            else:
                beta, p_wald = None, None
            
            data[f'beta_{trait_name}'] = beta
            data[f'p_wald_{trait_name}'] = p_wald
        
        rows.append(data)
    
    # Create DataFrame
    df = pd.DataFrame(rows)
    
    # Order columns
    base_cols = ['rs', 'chr', 'pos', 'allele1', 'allele0', 'af', 'trait_significant']
    effect_cols = []
    for trait_name in trait_names:
        effect_cols.extend([f'beta_{trait_name}', f'p_wald_{trait_name}'])
    
    df = df[base_cols + effect_cols]
    
    if verbose:
        print(f"  Created cross-trait table with {len(df):,} SNPs")
    
    return df

def main():
    """Main execution function."""
    args = parse_arguments()
    
    # Validate inputs
    try:
        validate_inputs(args.gwas_tables, args.trait_names, args.significant_snps)
    except Exception as e:
        print(f"Input validation error: {str(e)}")
        sys.exit(1)
    
    if args.verbose:
        print(f"Cross-Trait Effect Extraction")
        print(f"============================")
        print(f"GWAS traits: {', '.join(args.trait_names)}")
        print(f"Significant SNPs file: {args.significant_snps}")
        print(f"Output file: {args.output}")
        print()
    
    # Load GWAS data for all traits
    gwas_dicts = []
    for filepath, trait_name in zip(args.gwas_tables, args.trait_names):
        try:
            gwas_dict = load_gwas_data(filepath, trait_name, verbose=args.verbose)
            gwas_dicts.append(gwas_dict)
        except Exception as e:
            print(f"Error loading GWAS data for {trait_name} ({filepath}): {str(e)}")
            sys.exit(1)
    
    # Load significant SNPs
    try:
        sig_snps_df = load_significant_snps(args.significant_snps, verbose=args.verbose)
    except Exception as e:
        print(f"Error loading significant SNPs: {str(e)}")
        sys.exit(1)
    
    # Create cross-trait table
    try:
        cross_df = create_cross_trait_table(
            sig_snps_df=sig_snps_df,
            gwas_dicts=gwas_dicts,
            trait_names=args.trait_names,
            verbose=args.verbose
        )
    except Exception as e:
        print(f"Error creating cross-trait table: {str(e)}")
        sys.exit(1)
    
    # Save output
    try:
        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        cross_df.to_csv(output_path, sep='\t', index=False)
        
        if args.verbose:
            print(f"\nSaved cross-trait effects to: {output_path}")
    except Exception as e:
        print(f"Error saving output file: {str(e)}")
        sys.exit(1)
    
    # Generate summary statistics
    summary_data = []
    
    # Count SNPs per significance pattern
    for traits in cross_df['trait_significant'].unique():
        count = len(cross_df[cross_df['trait_significant'] == traits])
        summary_data.append({
            'significance_pattern': traits,
            'count': count,
            'percentage': 100 * count / len(cross_df)
        })
    
    # Create summary DataFrame
    summary_df = pd.DataFrame(summary_data)
    summary_df = summary_df.sort_values('count', ascending=False)
    
    # Save summary
    summary_path = Path(args.output).with_suffix('.summary.tsv')
    summary_df.to_csv(summary_path, sep='\t', index=False)
    
    if args.verbose:
        print(f"\nSummary statistics:")
        print("-" * 50)
        for _, row in summary_df.iterrows():
            print(f"{row['significance_pattern']:20s}: {row['count']:4d} SNPs ({row['percentage']:.1f}%)")
        print("-" * 50)
        print(f"Total SNPs in table: {len(cross_df)}")
        print(f"Summary saved to: {summary_path}")
    
    # Additional analysis: effect direction concordance
    if len(args.trait_names) >= 2:
        concordance_data = []
        
        # For each SNP, check effect direction concordance between traits
        for i, trait1 in enumerate(args.trait_names):
            for j, trait2 in enumerate(args.trait_names[i+1:], i+1):
                beta_col1 = f'beta_{trait1}'
                beta_col2 = f'beta_{trait2}'
                
                # Filter SNPs with effects in both traits
                valid_mask = cross_df[beta_col1].notna() & cross_df[beta_col2].notna()
                valid_df = cross_df[valid_mask].copy()
                
                if len(valid_df) > 0:
                    # Count concordant vs discordant
                    concordant = ((valid_df[beta_col1] > 0) & (valid_df[beta_col2] > 0)) | \
                                 ((valid_df[beta_col1] < 0) & (valid_df[beta_col2] < 0))
                    
                    n_concordant = concordant.sum()
                    n_discordant = len(valid_df) - n_concordant
                    
                    concordance_data.append({
                        'trait_pair': f'{trait1}_vs_{trait2}',
                        'n_snps': len(valid_df),
                        'n_concordant': n_concordant,
                        'n_discordant': n_discordant,
                        'concordance_rate': n_concordant / len(valid_df) if len(valid_df) > 0 else 0
                    })
        
        if concordance_data:
            concordance_df = pd.DataFrame(concordance_data)
            concordance_path = Path(args.output).with_suffix('.concordance.tsv')
            concordance_df.to_csv(concordance_path, sep='\t', index=False)
            
            if args.verbose:
                print(f"\nEffect direction concordance:")
                print("-" * 50)
                for _, row in concordance_df.iterrows():
                    print(f"{row['trait_pair']:20s}: {row['concordance_rate']:.1%} "
                          f"({row['n_concordant']}/{row['n_snps']} SNPs)")
                print("-" * 50)
                print(f"Concordance saved to: {concordance_path}")
    
    if args.verbose:
        print(f"\nAnalysis completed successfully")

if __name__ == "__main__":
    main()