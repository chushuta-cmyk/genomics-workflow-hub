#!/usr/bin/env python3
"""
Extract genome-wide significant SNPs from GWAS results.

This script filters SNPs with p_wald < threshold (default: 5e-8) from multiple
GWAS trait files and creates both trait-specific and merged significant SNP files.

Usage:
    python extract_significant_snps.py --gwas-tables size.tsv protein.tsv oil.tsv --trait-names size protein oil --output-dir SNP_analyze

Input:
    Simplified GWAS TSV files from extract_gwas_tables.py

Output:
    sig_<trait>.tsv - Significant SNPs for each trait
    sig_all_traits.tsv - All significant SNPs with trait annotation
"""

import pandas as pd
import argparse
import os
import sys
from pathlib import Path

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Extract genome-wide significant SNPs from GWAS results"
    )
    parser.add_argument(
        "--gwas-tables",
        nargs="+",
        required=True,
        help="List of GWAS table files (TSV format)"
    )
    parser.add_argument(
        "--trait-names",
        nargs="+",
        required=True,
        help="Corresponding trait names for each GWAS table"
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory for output files"
    )
    parser.add_argument(
        "--p-threshold",
        type=float,
        default=5e-8,
        help="P-value threshold for significance (default: 5e-8)"
    )
    parser.add_argument(
        "--min-af",
        type=float,
        default=0.01,
        help="Minimum allele frequency filter (default: 0.01)"
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print detailed progress information"
    )
    
    return parser.parse_args()

def validate_inputs(gwas_files, trait_names):
    """Validate input files and parameters."""
    if len(gwas_files) != len(trait_names):
        raise ValueError(f"Number of GWAS files ({len(gwas_files)}) does not match number of trait names ({len(trait_names)})")
    
    for filepath in gwas_files:
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"GWAS file not found: {filepath}")
    
    # Check for duplicate trait names
    if len(trait_names) != len(set(trait_names)):
        duplicates = [trait for trait in trait_names if trait_names.count(trait) > 1]
        raise ValueError(f"Duplicate trait names: {duplicates}")

def load_gwas_table(filepath, trait_name, p_threshold, min_af, verbose=False):
    """
    Load GWAS table and filter significant SNPs.
    
    Args:
        filepath: Path to GWAS TSV file
        trait_name: Name of the trait
        p_threshold: P-value threshold
        min_af: Minimum allele frequency
        verbose: Print progress if True
    
    Returns:
        DataFrame with significant SNPs
    """
    if verbose:
        print(f"Loading {trait_name} from {filepath}")
    
    # Load GWAS table
    df = pd.read_csv(filepath, sep='\t')
    
    # Check required columns
    required_cols = ['rs', 'chr', 'ps', 'allele1', 'allele0', 'af', 'beta', 'p_wald']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns in {filepath}: {missing_cols}")
    
    # Filter by p-value and allele frequency
    n_before = len(df)
    sig_df = df[(df['p_wald'] < p_threshold) & (df['af'] >= min_af) & (df['af'] <= 1 - min_af)].copy()
    n_after = len(sig_df)
    
    # Add trait column
    sig_df['trait'] = trait_name
    
    if verbose:
        print(f"  Total SNPs: {n_before:,}")
        print(f"  Significant SNPs (p < {p_threshold:.1e}): {n_after:,}")
        print(f"  Significant SNPs with AF >= {min_af}: {n_after:,}")
    
    return sig_df

def main():
    """Main execution function."""
    args = parse_arguments()
    
    # Validate inputs
    try:
        validate_inputs(args.gwas_tables, args.trait_names)
    except Exception as e:
        print(f"Input validation error: {str(e)}")
        sys.exit(1)
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    if args.verbose:
        print(f"Significant SNP Extraction")
        print(f"========================")
        print(f"P-value threshold: {args.p_threshold:.1e}")
        print(f"Minimum AF: {args.min_af}")
        print(f"Output directory: {output_dir}")
        print(f"GWAS files: {len(args.gwas_tables)}")
        print()
    
    # Process each trait
    all_sig_snps = []
    trait_summaries = []
    
    for filepath, trait_name in zip(args.gwas_tables, args.trait_names):
        try:
            # Load and filter significant SNPs
            sig_df = load_gwas_table(
                filepath=filepath,
                trait_name=trait_name,
                p_threshold=args.p_threshold,
                min_af=args.min_af,
                verbose=args.verbose
            )
            
            # Save trait-specific significant SNPs
            trait_output = output_dir / f"sig_{trait_name}.tsv"
            sig_df.to_csv(trait_output, sep='\t', index=False)
            
            if args.verbose:
                print(f"  Saved to: {trait_output}")
            
            # Add to combined list
            all_sig_snps.append(sig_df)
            
            # Store summary
            trait_summaries.append({
                'trait': trait_name,
                'total_snps': len(sig_df),
                'file': trait_output
            })
            
        except Exception as e:
            print(f"Error processing {trait_name} ({filepath}): {str(e)}")
            sys.exit(1)
    
    # Combine all significant SNPs
    if all_sig_snps:
        combined_df = pd.concat(all_sig_snps, ignore_index=True)
        
        # Save combined file
        combined_output = output_dir / "sig_all_traits.tsv"
        combined_df.to_csv(combined_output, sep='\t', index=False)
        
        if args.verbose:
            print(f"\nCombined significant SNPs: {len(combined_df):,}")
            print(f"Saved to: {combined_output}")
    
    # Generate summary statistics
    summary_output = output_dir / "sig_summary.txt"
    with open(summary_output, 'w') as f:
        f.write("Significant SNP Summary\n")
        f.write("======================\n")
        f.write(f"P-value threshold: {args.p_threshold:.1e}\n")
        f.write(f"Minimum AF: {args.min_af}\n")
        f.write(f"Analysis date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        f.write("Per-trait counts:\n")
        f.write("-" * 40 + "\n")
        for summary in trait_summaries:
            f.write(f"{summary['trait']:15s}: {summary['total_snps']:6d} SNPs\n")
        
        f.write("\n" + "-" * 40 + "\n")
        f.write(f"Total unique SNPs: {len(combined_df):6d}\n")
        
        # Count unique SNPs across traits
        if not combined_df.empty:
            unique_snps = combined_df['rs'].nunique()
            f.write(f"Unique SNP IDs:    {unique_snps:6d}\n")
            
            # SNPs significant in multiple traits
            snp_counts = combined_df['rs'].value_counts()
            multi_trait_snps = sum(snp_counts > 1)
            f.write(f"SNPs in >1 trait:  {multi_trait_snps:6d}\n")
    
    if args.verbose:
        print(f"\nSummary saved to: {summary_output}")
        
        # Print summary to console
        print("\n" + "=" * 50)
        print("SIGNIFICANT SNP SUMMARY")
        print("=" * 50)
        for summary in trait_summaries:
            print(f"{summary['trait']:15s}: {summary['total_snps']:6d} SNPs")
        print("-" * 50)
        
        if not combined_df.empty:
            unique_snps = combined_df['rs'].nunique()
            snp_counts = combined_df['rs'].value_counts()
            multi_trait_snps = sum(snp_counts > 1)
            
            print(f"Total unique SNPs: {len(combined_df):6d}")
            print(f"Unique SNP IDs:    {unique_snps:6d}")
            print(f"SNPs in >1 trait:  {multi_trait_snps:6d}")
        
        print("=" * 50)
        print(f"\nAnalysis completed successfully")
        print(f"Output directory: {output_dir}")

if __name__ == "__main__":
    main()