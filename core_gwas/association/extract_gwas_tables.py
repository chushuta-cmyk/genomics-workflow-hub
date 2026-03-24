#!/usr/bin/env python3
"""
Extract essential columns from GWAS result files.

This script reads GWAS association files (GEMMA .assoc.txt format) and extracts
essential columns for downstream analysis. It creates simplified TSV files
with only the necessary columns for efficient processing.

Usage:
    python extract_gwas_tables.py --index gwas_files.txt --output-dir SNP_analyze

Input:
    GWAS association files in GEMMA format with columns:
    chr, rs, ps, n_miss, allele1, allele0, af, beta, se, logl_H1, l_remle, l_mle, p_wald, p_lrt, p_score

Output:
    Simplified TSV files with columns:
    chr, rs, ps, allele1, allele0, af, beta, se, p_wald
"""

import pandas as pd
import argparse
import os
import sys
from pathlib import Path

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Extract essential columns from GWAS result files"
    )
    parser.add_argument(
        "--index",
        required=True,
        help="Path to GWAS file index (format: trait:path per line)"
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory for output files"
    )
    parser.add_argument(
        "--columns",
        nargs="+",
        default=["chr", "rs", "ps", "allele1", "allele0", "af", "beta", "se", "p_wald"],
        help="Columns to extract (default: chr rs ps allele1 allele0 af beta se p_wald)"
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print detailed progress information"
    )
    parser.add_argument(
        "--chunksize",
        type=int,
        default=None,
        help="Process large files in chunks (number of rows per chunk)"
    )
    
    return parser.parse_args()

def read_gwas_index(index_file):
    """
    Read GWAS file index.
    
    Format: trait_name:/path/to/gwas_file.assoc.txt
    One per line, trait name and path separated by colon.
    
    Returns:
        Dictionary {trait_name: file_path}
    """
    gwas_files = {}
    with open(index_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                parts = line.split(':', 1)
                if len(parts) == 2:
                    trait = parts[0].strip()
                    path = parts[1].strip()
                    gwas_files[trait] = path
                else:
                    print(f"Warning: Skipping malformed line: {line}")
    
    return gwas_files

def extract_gwas_columns(input_file, output_file, columns, chunksize=None, verbose=False):
    """
    Extract specified columns from GWAS file.
    
    Args:
        input_file: Path to GWAS association file
        output_file: Path for output TSV file
        columns: List of column names to extract
        chunksize: Process in chunks if specified
        verbose: Print progress if True
    """
    if verbose:
        print(f"Processing: {input_file}")
        print(f"  Output: {output_file}")
        print(f"  Columns: {', '.join(columns)}")
    
    # Check if input file exists
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file not found: {input_file}")
    
    # Read GWAS file
    try:
        if chunksize:
            # Process in chunks for large files
            chunks = []
            for i, chunk in enumerate(pd.read_csv(input_file, sep='\t', chunksize=chunksize)):
                # Select required columns
                chunk = chunk[columns]
                chunks.append(chunk)
                if verbose and i % 10 == 0:
                    print(f"  Processed {i * chunksize:,} rows...")
            
            df = pd.concat(chunks, ignore_index=True)
        else:
            # Read entire file
            df = pd.read_csv(input_file, sep='\t')
            # Select required columns
            df = df[columns]
        
        # Ensure rs column is string
        df['rs'] = df['rs'].astype(str)
        
        # Save to output file
        df.to_csv(output_file, sep='\t', index=False)
        
        if verbose:
            print(f"  Saved {len(df):,} SNPs to {output_file}")
        
        return len(df)
        
    except Exception as e:
        print(f"Error processing {input_file}: {str(e)}")
        raise

def main():
    """Main execution function."""
    args = parse_arguments()
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Read GWAS file index
    try:
        gwas_files = read_gwas_index(args.index)
        if not gwas_files:
            print(f"Error: No valid GWAS files found in index: {args.index}")
            sys.exit(1)
        
        if args.verbose:
            print(f"Found {len(gwas_files)} GWAS files in index:")
            for trait, path in gwas_files.items():
                print(f"  {trait}: {path}")
    
    except Exception as e:
        print(f"Error reading index file {args.index}: {str(e)}")
        sys.exit(1)
    
    # Process each GWAS file
    total_snps = 0
    success_count = 0
    
    for trait, input_file in gwas_files.items():
        try:
            # Create output filename
            output_file = output_dir / f"gwas_{trait}.tsv"
            
            # Extract columns
            n_snps = extract_gwas_columns(
                input_file=input_file,
                output_file=output_file,
                columns=args.columns,
                chunksize=args.chunksize,
                verbose=args.verbose
            )
            
            total_snps += n_snps
            success_count += 1
            
        except Exception as e:
            print(f"Failed to process {trait} ({input_file}): {str(e)}")
    
    # Summary
    if args.verbose:
        print(f"\nSummary:")
        print(f"  Successfully processed: {success_count}/{len(gwas_files)} files")
        print(f"  Total SNPs extracted: {total_snps:,}")
        print(f"  Output directory: {output_dir}")
    
    # Check if all files were processed
    if success_count == 0:
        print("Error: No files were successfully processed")
        sys.exit(1)
    elif success_count < len(gwas_files):
        print(f"Warning: Only {success_count}/{len(gwas_files)} files processed successfully")
        sys.exit(0)
    else:
        if args.verbose:
            print("All files processed successfully")
        sys.exit(0)

if __name__ == "__main__":
    main()