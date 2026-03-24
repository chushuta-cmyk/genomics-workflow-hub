#!/usr/bin/env python3
"""
Generate final analysis summary report.

This script consolidates results from all analysis steps into a comprehensive
markdown report with tables, statistics, and interpretation.

Usage:
    python generate_summary.py --sig-snps sig_all_traits.tsv --cross-trait cross_trait_effect.tsv --loci locus_summary.tsv --coloc coloc_results.tsv --output analysis_summary.md

Input:
    Results from all previous analysis steps

Output:
    analysis_summary.md - Comprehensive markdown report
"""

import pandas as pd
import argparse
import sys
from pathlib import Path
from datetime import datetime

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate final analysis summary report"
    )
    parser.add_argument(
        "--sig-snps",
        required=True,
        help="Significant SNPs file (sig_all_traits.tsv)"
    )
    parser.add_argument(
        "--cross-trait",
        required=True,
        help="Cross-trait effects file (cross_trait_effect.tsv)"
    )
    parser.add_argument(
        "--loci",
        required=True,
        help="Locus summary file (locus_summary.tsv)"
    )
    parser.add_argument(
        "--coloc",
        required=True,
        help="Colocalization results file (coloc_results.tsv)"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output markdown file path"
    )
    parser.add_argument(
        "--title",
        default="GWAS Post-Analysis Summary Report",
        help="Report title (default: GWAS Post-Analysis Summary Report)"
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print detailed progress information"
    )
    
    return parser.parse_args()

def load_data(sig_snps_file, cross_trait_file, loci_file, coloc_file, verbose=False):
    """Load all input data files."""
    data = {}
    
    try:
        if verbose:
            print("Loading data files...")
        
        # Load significant SNPs
        data['sig_snps'] = pd.read_csv(sig_snps_file, sep='\t')
        if verbose:
            print(f"  Significant SNPs: {len(data['sig_snps']):,} rows")
        
        # Load cross-trait effects
        data['cross_trait'] = pd.read_csv(cross_trait_file, sep='\t')
        if verbose:
            print(f"  Cross-trait effects: {len(data['cross_trait']):,} rows")
        
        # Load locus summary
        data['loci'] = pd.read_csv(loci_file, sep='\t')
        if verbose:
            print(f"  Loci: {len(data['loci']):,} rows")
        
        # Load colocalization results
        data['coloc'] = pd.read_csv(coloc_file, sep='\t')
        if verbose:
            print(f"  Colocalization results: {len(data['coloc']):,} rows")
        
    except Exception as e:
        print(f"Error loading data: {str(e)}")
        sys.exit(1)
    
    return data

def generate_summary_statistics(data, verbose=False):
    """Generate comprehensive summary statistics."""
    stats = {}
    
    # 1. Significant SNP statistics
    sig_snps = data['sig_snps']
    if 'trait' in sig_snps.columns:
        stats['sig_by_trait'] = sig_snps['trait'].value_counts().to_dict()
        stats['total_sig_snps'] = len(sig_snps)
        stats['unique_sig_snps'] = sig_snps['rs'].nunique()
        
        # SNPs significant in multiple traits
        snp_counts = sig_snps['rs'].value_counts()
        stats['multi_trait_snps'] = sum(snp_counts > 1)
        stats['single_trait_snps'] = sum(snp_counts == 1)
    
    # 2. Cross-trait statistics
    cross_trait = data['cross_trait']
    if 'trait_significant' in cross_trait.columns:
        trait_patterns = cross_trait['trait_significant'].value_counts()
        stats['trait_patterns'] = trait_patterns.head(10).to_dict()  # Top 10 patterns
    
    # 3. Locus statistics
    loci = data['loci']
    if not loci.empty:
        stats['total_loci'] = len(loci)
        
        if 'num_snps' in loci.columns:
            stats['loci_with_multiple_snps'] = len(loci[loci['num_snps'] > 1])
            stats['largest_locus_snps'] = loci['num_snps'].max()
        
        if 'locus_size' in loci.columns:
            stats['largest_locus_size'] = loci['locus_size'].max()
        
        # Loci by trait
        if 'traits_in_locus' in loci.columns:
            multi_trait_loci = loci[loci['traits_in_locus'].str.contains(',', na=False)]
            stats['multi_trait_loci'] = len(multi_trait_loci)
            stats['single_trait_loci'] = len(loci) - len(multi_trait_loci)
    
    # 4. Colocalization statistics
    coloc = data['coloc']
    if not coloc.empty and 'coloc_colocalized' in coloc.columns:
        stats['total_coloc_tests'] = len(coloc)
        stats['colocalized_loci'] = coloc['coloc_colocalized'].sum()
        stats['colocalization_rate'] = coloc['coloc_colocalized'].mean()
        
        # By trait pair
        if 'trait_pair' in coloc.columns:
            coloc_by_pair = coloc.groupby('trait_pair')['coloc_colocalized'].agg(['count', 'sum', 'mean'])
            stats['coloc_by_pair'] = coloc_by_pair.to_dict('index')
    
    return stats

def create_markdown_report(data, stats, args):
    """Create markdown report with all sections."""
    report = []
    
    # Header
    report.append(f"# {args.title}")
    report.append("")
    report.append(f"**Date:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    report.append(f"**Generated by:** GWAS Post-Analysis Pipeline v1.0")
    report.append("")
    report.append("---")
    report.append("")
    
    # Executive Summary
    report.append("## Executive Summary")
    report.append("")
    report.append("This report summarizes the results of the multi-trait GWAS post-analysis pipeline, including significant SNP extraction, cross-trait effect comparison, locus merging, and Bayesian colocalization analysis.")
    report.append("")
    
    # Key Findings
    report.append("### Key Findings")
    report.append("")
    
    if 'total_sig_snps' in stats:
        report.append(f"1. **{stats['total_sig_snps']:,} genome-wide significant SNPs** identified across all traits (p < 5e-8)")
    
    if 'unique_sig_snps' in stats:
        report.append(f"2. **{stats['unique_sig_snps']:,} unique SNP IDs** with genome-wide significance")
    
    if 'multi_trait_snps' in stats:
        report.append(f"3. **{stats['multi_trait_snps']:,} SNPs** are significant for multiple traits")
    
    if 'total_loci' in stats:
        report.append(f"4. **{stats['total_loci']:,} genomic loci** identified (±250kb window)")
    
    if 'colocalized_loci' in stats and 'total_coloc_tests' in stats:
        report.append(f"5. **{stats['colocalized_loci']:,} of {stats['total_coloc_tests']:,} loci** show evidence of colocalization (PP4 ≥ 0.8)")
    
    report.append("")
    report.append("---")
    report.append("")
    
    # 1. Significant SNP Analysis
    report.append("## 1. Significant SNP Analysis")
    report.append("")
    report.append("### 1.1 Per-Trait Counts")
    report.append("")
    
    if 'sig_by_trait' in stats:
        report.append("| Trait | Significant SNPs |")
        report.append("|-------|------------------|")
        for trait, count in stats['sig_by_trait'].items():
            report.append(f"| {trait} | {count:,} |")
        report.append("")
        
        report.append(f"**Total:** {stats['total_sig_snps']:,} significant SNP-trait associations")
        report.append("")
        report.append(f"**Unique SNPs:** {stats['unique_sig_snps']:,}")
        report.append("")
        report.append(f"**SNPs in multiple traits:** {stats['multi_trait_snps']:,}")
        report.append("")
    
    # 2. Cross-Trait Effects
    report.append("## 2. Cross-Trait Effects")
    report.append("")
    report.append("### 2.1 Trait Significance Patterns")
    report.append("")
    
    if 'trait_patterns' in stats:
        report.append("| Trait Pattern | Count | Percentage |")
        report.append("|---------------|-------|------------|")
        total_patterns = sum(stats['trait_patterns'].values())
        for pattern, count in stats['trait_patterns'].items():
            percentage = 100 * count / total_patterns
            report.append(f"| {pattern} | {count:,} | {percentage:.1f}% |")
        report.append("")
    
    # 3. Locus Analysis
    report.append("## 3. Locus Analysis")
    report.append("")
    report.append("### 3.1 Locus Summary Statistics")
    report.append("")
    
    if 'total_loci' in stats:
        report.append(f"- **Total loci identified:** {stats['total_loci']:,}")
        report.append(f"- **Loci with multiple SNPs:** {stats['loci_with_multiple_snps']:,}")
        
        if 'multi_trait_loci' in stats:
            report.append(f"- **Multi-trait loci:** {stats['multi_trait_loci']:,}")
            report.append(f"- **Single-trait loci:** {stats['single_trait_loci']:,}")
        
        if 'largest_locus_snps' in stats:
            report.append(f"- **Largest locus (SNP count):** {stats['largest_locus_snps']:,} SNPs")
        
        if 'largest_locus_size' in stats:
            report.append(f"- **Largest locus (genomic span):** {stats['largest_locus_size']:,} bp")
        report.append("")
    
    # 4. Colocalization Analysis
    report.append("## 4. Colocalization Analysis")
    report.append("")
    report.append("### 4.1 Overall Results")
    report.append("")
    
    if 'total_coloc_tests' in stats and 'colocalized_loci' in stats:
        report.append(f"- **Total colocalization tests:** {stats['total_coloc_tests']:,}")
        report.append(f"- **Colocalized loci:** {stats['colocalized_loci']:,}")
        report.append(f"- **Colocalization rate:** {100 * stats['colocalization_rate']:.1f}%")
        report.append("")
        
        report.append("*Interpretation:* PP4 ≥ 0.8 indicates strong evidence for shared causal variants between trait pairs at the same genomic locus.")
        report.append("")
    
    # 4.2 By Trait Pair
    if 'coloc_by_pair' in stats:
        report.append("### 4.2 Colocalization by Trait Pair")
        report.append("")
        report.append("| Trait Pair | Tests | Colocalized | Rate |")
        report.append("|------------|-------|-------------|------|")
        
        for trait_pair, pair_stats in stats['coloc_by_pair'].items():
            rate = pair_stats['mean']
            report.append(f"| {trait_pair} | {pair_stats['count']:,} | {pair_stats['sum']:,} | {100 * rate:.1f}% |")
        report.append("")
    
    # 5. Top Loci
    report.append("## 5. Top Loci of Interest")
    report.append("")
    
    # Identify top loci by various criteria
    loci = data['loci']
    if not loci.empty:
        # Top loci by number of SNPs
        if 'num_snps' in loci.columns:
            top_snp_loci = loci.nlargest(5, 'num_snps')[['locus_id', 'chr', 'start', 'end', 'num_snps', 'traits_in_locus']]
            
            report.append("### 5.1 Loci with Most SNPs")
            report.append("")
            report.append("| Locus ID | Chromosome | Start | End | SNPs | Traits |")
            report.append("|----------|------------|-------|-----|------|--------|")
            for _, row in top_snp_loci.iterrows():
                report.append(f"| {row['locus_id']} | {row['chr']} | {row['start']:,} | {row['end']:,} | {row['num_snps']} | {row['traits_in_locus']} |")
            report.append("")
        
        # Top loci by colocalization evidence
        coloc = data['coloc']
        if not coloc.empty and 'coloc_pp4' in coloc.columns:
            top_coloc_loci = coloc.nlargest(5, 'coloc_pp4')[['locus_id', 'trait_pair', 'coloc_pp4', 'coloc_colocalized']]
            
            report.append("### 5.2 Loci with Strongest Colocalization Evidence")
            report.append("")
            report.append("| Locus ID | Trait Pair | PP4 | Colocalized |")
            report.append("|----------|------------|-----|-------------|")
            for _, row in top_coloc_loci.iterrows():
                report.append(f"| {row['locus_id']} | {row['trait_pair']} | {row['coloc_pp4']:.3f} | {'Yes' if row['coloc_colocalized'] else 'No'} |")
            report.append("")
    
    # 6. Methods
    report.append("## 6. Methods")
    report.append("")
    report.append("### 6.1 Analysis Pipeline")
    report.append("")
    report.append("The analysis followed these steps:")
    report.append("")
    report.append("1. **GWAS Result Indexing**: Identified and indexed GWAS result files")
    report.append("2. **GWAS Table Extraction**: Extracted essential columns from GWAS files")
    report.append("3. **Significant SNP Extraction**: Filtered SNPs with p_wald < 5e-8")
    report.append("4. **Cross-Trait Effect Extraction**: Retrieved SNP effects across all traits")
    report.append("5. **Locus Merging**: Grouped SNPs within ±250kb into loci")
    report.append("6. **Colocalization Analysis**: Applied Bayesian colocalization (coloc) to test for shared causal variants")
    report.append("7. **Report Generation**: Compiled results into this comprehensive summary")
    report.append("")
    
    report.append("### 6.2 Statistical Thresholds")
    report.append("")
    report.append("- **Genome-wide significance**: p < 5e-8")
    report.append("- **Locus definition**: SNPs within ±250kb merged into same locus")
    report.append("- **Colocalization threshold**: PP4 ≥ 0.8")
    report.append("- **Allele frequency filter**: MAF ≥ 0.01")
    report.append("")
    
    # 7. Files Generated
    report.append("## 7. Files Generated")
    report.append("")
    report.append("| File | Description |")
    report.append("|------|-------------|")
    
    # List input files
    report.append(f"| `{args.sig_snps}` | Significant SNPs across all traits |")
    report.append(f"| `{args.cross_trait}` | Cross-trait SNP effects |")
    report.append(f"| `{args.loci}` | Merged locus summary |")
    report.append(f"| `{args.coloc}` | Colocalization results |")
    
    # Add derived files
    base_output = Path(args.output)
    report.append(f"| `{base_output}` | This summary report |")
    
    # Common derived files
    derived_files = [
        "sig_summary.txt",
        "cross_trait_effect.summary.tsv",
        "cross_trait_effect.concordance.tsv",
        "locus_summary.stats.tsv",
        "coloc_results.summary.tsv"
    ]
    
    for file in derived_files:
        derived_path = base_output.parent / file
        if derived_path.exists():
            report.append(f"| `{file}` | Derived statistics |")
    
    report.append("")
    
    # 8. Limitations
    report.append("## 8. Limitations and Considerations")
    report.append("")
    report.append("1. **Sample size**: GWAS power depends on sample size and phenotypic variance")
    report.append("2. **Population stratification**: Results may be population-specific")
    report.append("3. **LD structure**: Locus merging depends on local LD patterns")
    report.append("4. **Colocalization assumptions**: Methods assume single causal variant per trait per locus")
    report.append("5. **Multiple testing**: Consideration needed for the number of loci and trait pairs tested")
    report.append("")
    
    # 9. Conclusion
    report.append("## 9. Conclusion")
    report.append("")
    report.append("This analysis provides a comprehensive overview of genetic associations across multiple traits, identifying:")
    report.append("")
    report.append("1. **Trait-specific and pleiotropic SNPs** contributing to phenotypic variation")
    report.append("2. **Genomic loci** harboring multiple associated variants")
    report.append("3. **Evidence for shared genetic architecture** through colocalization analysis")
    report.append("")
    report.append("These results can inform follow-up studies including fine-mapping, functional validation, and breeding applications.")
    report.append("")
    
    # Footer
    report.append("---")
    report.append("")
    report.append("*Report generated automatically by GWAS Post-Analysis Pipeline*")
    report.append("*For questions or further analysis, contact the bioinformatics team*")
    
    return "\n".join(report)

def main():
    """Main execution function."""
    args = parse_arguments()
    
    if args.verbose:
        print(f"Generating Analysis Summary Report")
        print(f"==================================")
        print(f"Title: {args.title}")
        print(f"Output: {args.output}")
        print()
    
    # Load data
    data = load_data(
        sig_snps_file=args.sig_snps,
        cross_trait_file=args.cross_trait,
        loci_file=args.loci,
        coloc_file=args.coloc,
        verbose=args.verbose
    )
    
    # Generate statistics
    stats = generate_summary_statistics(data, verbose=args.verbose)
    
    # Create markdown report
    report = create_markdown_report(data, stats, args)
    
    # Save report
    try:
        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_path, 'w') as f:
            f.write(report)
        
        if args.verbose:
            print(f"Report saved to: {output_path}")
            print(f"Report length: {len(report):,} characters")
            print(f"Report sections: {report.count('## ')} main sections")
            print()
            print("Analysis summary completed successfully")
    
    except Exception as e:
        print(f"Error saving report: {str(e)}")
        sys.exit(1)
    
    # Also save statistics as JSON for programmatic access
    if args.verbose:
        import json
        stats_path = output_path.with_suffix('.stats.json')
        with open(stats_path, 'w') as f:
            json.dump(stats, f, indent=2, default=str)
        print(f"Statistics saved to: {stats_path}")

if __name__ == "__main__":
    main()