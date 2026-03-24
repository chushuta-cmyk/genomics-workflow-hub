#!/usr/bin/env python3
"""
Main pipeline runner for GWAS post-analysis.

This script orchestrates the complete 7-step workflow based on configuration.

Usage:
    python run_pipeline.py --config config_analysis.yaml

Or for quick execution:
    python run_pipeline.py --gwas-files size.assoc.txt protein.assoc.txt oil.assoc.txt --trait-names size protein oil
"""

import argparse
import os
import sys
import yaml
import subprocess
from pathlib import Path
from datetime import datetime

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="GWAS Post-Analysis Pipeline Runner"
    )
    
    # Configuration file option
    parser.add_argument(
        "--config",
        help="Path to configuration file (YAML format)"
    )
    
    # Quick execution options (alternative to config file)
    parser.add_argument(
        "--gwas-files",
        nargs="+",
        help="List of GWAS association files"
    )
    parser.add_argument(
        "--trait-names",
        nargs="+",
        help="Corresponding trait names"
    )
    parser.add_argument(
        "--output-dir",
        default="./results",
        help="Output directory (default: ./results)"
    )
    
    # Pipeline control
    parser.add_argument(
        "--steps",
        nargs="+",
        default=["all"],
        help="Steps to execute: index, extract, significant, cross, loci, coloc, summary, all"
    )
    
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print detailed progress information"
    )
    
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands without executing"
    )
    
    return parser.parse_args()

def load_configuration(config_file):
    """Load configuration from YAML file."""
    try:
        with open(config_file, 'r') as f:
            config = yaml.safe_load(f)
        return config
    except Exception as e:
        print(f"Error loading configuration file {config_file}: {str(e)}")
        sys.exit(1)

def create_quick_config(args):
    """Create configuration dictionary from command line arguments."""
    config = {
        'dataset': {
            'type': 'custom',
            'name': 'quick_analysis'
        },
        'traits': {
            'names': args.trait_names
        },
        'analysis': {
            'p_threshold': 5e-8,
            'min_af': 0.01,
            'locus': {
                'window': 250000,
                'p_column': f"p_wald_{args.trait_names[0]}"
            },
            'coloc': {
                'pp4_threshold': 0.8,
                'priors': {
                    'p1': 1e-4,
                    'p2': 1e-4,
                    'p12': 1e-5
                }
            }
        },
        'output': {
            'directories': {
                'base': args.output_dir,
                'analysis_index': "analysis_index",
                'snp_analyze': "SNP_analyze",
                'logs': "logs",
                'scripts': "scripts"
            }
        },
        'pipeline': {
            'steps': {
                'index_gwas_files': True,
                'extract_gwas_tables': True,
                'extract_significant_snps': True,
                'cross_trait_effects': True,
                'merge_loci': True,
                'run_coloc': True,
                'generate_summary': True
            }
        }
    }
    
    # Add GWAS file paths
    config['dataset']['gwas_files'] = {}
    for trait, filepath in zip(args.trait_names, args.gwas_files):
        config['dataset']['gwas_files'][trait] = filepath
    
    return config

def validate_configuration(config):
    """Validate configuration parameters."""
    required_sections = ['dataset', 'traits', 'analysis', 'output']
    for section in required_sections:
        if section not in config:
            raise ValueError(f"Missing required section: {section}")
    
    # Check GWAS files exist
    if 'gwas_files' in config['dataset']:
        for trait, filepath in config['dataset']['gwas_files'].items():
            if not os.path.exists(filepath):
                print(f"Warning: GWAS file for {trait} not found: {filepath}")
    
    return True

def setup_directories(config, verbose=False):
    """Create output directory structure."""
    base_dir = config['output']['directories']['base']
    directories = [
        base_dir,
        os.path.join(base_dir, config['output']['directories']['analysis_index']),
        os.path.join(base_dir, config['output']['directories']['snp_analyze']),
        os.path.join(base_dir, config['output']['directories']['logs']),
        os.path.join(base_dir, config['output']['directories']['scripts'])
    ]
    
    for directory in directories:
        os.makedirs(directory, exist_ok=True)
        if verbose:
            print(f"Created directory: {directory}")
    
    return base_dir

def run_command(cmd, log_file=None, verbose=False, dry_run=False):
    """Run a command and log output."""
    if dry_run:
        print(f"[DRY RUN] {cmd}")
        if log_file:
            print(f"  Log would be written to: {log_file}")
        return 0
    
    if verbose:
        print(f"Running: {cmd}")
    
    try:
        if log_file:
            with open(log_file, 'w') as log:
                result = subprocess.run(
                    cmd, shell=True, check=True,
                    stdout=log, stderr=subprocess.STDOUT
                )
            if verbose:
                print(f"  Output logged to: {log_file}")
        else:
            result = subprocess.run(cmd, shell=True, check=True, capture_output=True)
            if verbose and result.stdout:
                print(result.stdout.decode())
        
        return 0
    except subprocess.CalledProcessError as e:
        print(f"Command failed: {cmd}")
        print(f"Error code: {e.returncode}")
        if log_file and os.path.exists(log_file):
            with open(log_file, 'r') as f:
                print(f"Log contents:\n{f.read()}")
        return e.returncode

def step1_index_gwas_files(config, base_dir, verbose=False, dry_run=False):
    """Step 1: Create GWAS file index."""
    if verbose:
        print("\n=== Step 1: Indexing GWAS Files ===")
    
    index_file = os.path.join(
        base_dir,
        config['output']['directories']['analysis_index'],
        config['output']['files']['gwas_index']
    )
    
    gwas_files = config['dataset']['gwas_files']
    
    # Create index file
    index_content = []
    for trait, filepath in gwas_files.items():
        index_content.append(f"{trait}:{filepath}")
    
    if not dry_run:
        with open(index_file, 'w') as f:
            f.write('\n'.join(index_content))
    
    if verbose:
        print(f"Created index file: {index_file}")
        print(f"Indexed {len(gwas_files)} GWAS files")
    
    return index_file

def step2_extract_gwas_tables(config, base_dir, index_file, verbose=False, dry_run=False):
    """Step 2: Extract essential columns from GWAS files."""
    if verbose:
        print("\n=== Step 2: Extracting GWAS Tables ===")
    
    output_dir = os.path.join(base_dir, config['output']['directories']['snp_analyze'])
    log_file = os.path.join(base_dir, config['output']['directories']['logs'], "step2_extract.log")
    
    script_path = os.path.join(os.path.dirname(__file__), "scripts", "extract_gwas_tables.py")
    
    cmd = f"python {script_path} --index {index_file} --output-dir {output_dir} --verbose"
    
    return run_command(cmd, log_file, verbose, dry_run)

def step3_extract_significant_snps(config, base_dir, verbose=False, dry_run=False):
    """Step 3: Extract genome-wide significant SNPs."""
    if verbose:
        print("\n=== Step 3: Extracting Significant SNPs ===")
    
    output_dir = os.path.join(base_dir, config['output']['directories']['snp_analyze'])
    log_file = os.path.join(base_dir, config['output']['directories']['logs'], "step3_significant.log")
    
    # Get GWAS table files
    gwas_tables = []
    for trait in config['traits']['names']:
        table_file = os.path.join(output_dir, f"{config['output']['files']['extracted_prefix']}{trait}.tsv")
        gwas_tables.append(table_file)
    
    script_path = os.path.join(os.path.dirname(__file__), "scripts", "extract_significant_snps.py")
    
    cmd = f"python {script_path} --gwas-tables {' '.join(gwas_tables)} "
    cmd += f"--trait-names {' '.join(config['traits']['names'])} "
    cmd += f"--output-dir {output_dir} "
    cmd += f"--p-threshold {config['analysis']['p_threshold']} "
    cmd += f"--min-af {config['analysis']['min_af']} "
    cmd += "--verbose"
    
    return run_command(cmd, log_file, verbose, dry_run)

def step4_cross_trait_effects(config, base_dir, verbose=False, dry_run=False):
    """Step 4: Extract cross-trait effects."""
    if verbose:
        print("\n=== Step 4: Extracting Cross-Trait Effects ===")
    
    output_dir = os.path.join(base_dir, config['output']['directories']['snp_analyze'])
    log_file = os.path.join(base_dir, config['output']['directories']['logs'], "step4_cross_trait.log")
    
    # Get GWAS table files
    gwas_tables = []
    for trait in config['traits']['names']:
        table_file = os.path.join(output_dir, f"{config['output']['files']['extracted_prefix']}{trait}.tsv")
        gwas_tables.append(table_file)
    
    # Significant SNPs file
    sig_snps_file = os.path.join(output_dir, "sig_all_traits.tsv")
    
    # Output file
    cross_trait_file = os.path.join(output_dir, config['output']['files']['cross_trait'])
    
    script_path = os.path.join(os.path.dirname(__file__), "scripts", "cross_trait_effects.py")
    
    cmd = f"python {script_path} --gwas-tables {' '.join(gwas_tables)} "
    cmd += f"--trait-names {' '.join(config['traits']['names'])} "
    cmd += f"--significant-snps {sig_snps_file} "
    cmd += f"--output {cross_trait_file} "
    cmd += "--verbose"
    
    return run_command(cmd, log_file, verbose, dry_run)

def step5_merge_loci(config, base_dir, verbose=False, dry_run=False):
    """Step 5: Merge SNPs into loci."""
    if verbose:
        print("\n=== Step 5: Merging Loci ===")
    
    output_dir = os.path.join(base_dir, config['output']['directories']['snp_analyze'])
    log_file = os.path.join(base_dir, config['output']['directories']['logs'], "step5_loci.log")
    
    # Cross-trait effects file
    cross_trait_file = os.path.join(output_dir, config['output']['files']['cross_trait'])
    
    # Output file
    locus_file = os.path.join(output_dir, config['output']['files']['locus_summary'])
    
    script_path = os.path.join(os.path.dirname(__file__), "scripts", "merge_loci.py")
    
    cmd = f"python {script_path} --snp-table {cross_trait_file} "
    cmd += f"--window {config['analysis']['locus']['window']} "
    cmd += f"--p-column {config['analysis']['locus']['p_column']} "
    cmd += f"--output {locus_file} "
    cmd += "--verbose"
    
    return run_command(cmd, log_file, verbose, dry_run)

def step6_run_coloc(config, base_dir, verbose=False, dry_run=False):
    """Step 6: Run colocalization analysis."""
    if verbose:
        print("\n=== Step 6: Running Colocalization ===")
    
    output_dir = os.path.join(base_dir, config['output']['directories']['snp_analyze'])
    log_file = os.path.join(base_dir, config['output']['directories']['logs'], "step6_coloc.log")
    
    # Get GWAS table files
    gwas_tables = []
    for trait in config['traits']['names']:
        table_file = os.path.join(output_dir, f"{config['output']['files']['extracted_prefix']}{trait}.tsv")
        gwas_tables.append(table_file)
    
    # Locus file
    locus_file = os.path.join(output_dir, config['output']['files']['locus_summary'])
    
    # Output file
    coloc_file = os.path.join(output_dir, config['output']['files']['coloc_results'])
    
    script_path = os.path.join(os.path.dirname(__file__), "scripts", "run_coloc.R")
    
    cmd = f"Rscript {script_path} --gwas-files {' '.join(gwas_tables)} "
    cmd += f"--trait-names {' '.join(config['traits']['names'])} "
    cmd += f"--loci {locus_file} "
    cmd += f"--output {coloc_file} "
    cmd += f"--pp4-threshold {config['analysis']['coloc']['pp4_threshold']} "
    cmd += f"--prior-p1 {config['analysis']['coloc']['priors']['p1']} "
    cmd += f"--prior-p2 {config['analysis']['coloc']['priors']['p2']} "
    cmd += f"--prior-p12 {config['analysis']['coloc']['priors']['p12']} "
    cmd += "--verbose"
    
    return run_command(cmd, log_file, verbose, dry_run)

def step7_generate_summary(config, base_dir, verbose=False, dry_run=False):
    """Step 7: Generate final summary report."""
    if verbose:
        print("\n=== Step 7: Generating Summary Report ===")
    
    output_dir = os.path.join(base_dir, config['output']['directories']['snp_analyze'])
    log_file = os.path.join(base_dir, config['output']['directories']['logs'], "step7_summary.log")
    
    # Input files
    sig_snps_file = os.path.join(output_dir, "sig_all_traits.tsv")
    cross_trait_file = os.path.join(output_dir, config['output']['files']['cross_trait'])
    locus_file = os.path.join(output_dir, config['output']['files']['locus_summary'])
    coloc_file = os.path.join(output_dir, config['output']['files']['coloc_results'])
    
    # Output file
    summary_file = os.path.join(output_dir, config['output']['files']['analysis_summary'])
    
    script_path = os.path.join(os.path.dirname(__file__), "scripts", "generate_summary.py")
    
    cmd = f"python {script_path} --sig-snps {sig_snps_file} "
    cmd += f"--cross-trait {cross_trait_file} "
    cmd += f"--loci {locus_file} "
    cmd += f"--coloc {coloc_file} "
    cmd += f"--output {summary_file} "
    cmd += f"--title \"{config.get('metadata', {}).get('project', 'GWAS Post-Analysis Summary')}\" "
    cmd += "--verbose"
    
    return run_command(cmd, log_file, verbose, dry_run)

def main():
    """Main pipeline execution."""
    args = parse_arguments()
    
    # Load configuration
    if args.config:
        config = load_configuration(args.config)
    elif args.gwas_files and args.trait_names:
        if len(args.gwas_files) != len(args.trait_names):
            print("Error: Number of GWAS files must match number of trait names")
            sys.exit(1)
        config = create_quick_config(args)
    else:
        print("Error: Either --config or both --gwas-files and --trait-names must be provided")
        sys.exit(1)
    
    # Validate configuration
    try:
        validate_configuration(config)
    except Exception as e:
        print(f"Configuration error: {str(e)}")
        sys.exit(1)
    
    # Setup directories
    base_dir = setup_directories(config, args.verbose)
    
    # Determine which steps to run
    if "all" in args.steps:
        steps_to_run = [1, 2, 3, 4, 5, 6, 7]
    else:
        steps_to_run = []
        step_map = {
            "index": 1,
            "extract": 2,
            "significant": 3,
            "cross": 4,
            "loci": 5,
            "coloc": 6,
            "summary": 7
        }
        for step_name in args.steps:
            if step_name in step_map:
                steps_to_run.append(step_map[step_name])
    
    if args.verbose:
        print(f"Pipeline starting at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"Base directory: {base_dir}")
        print(f"Steps to execute: {steps_to_run}")
        print()
    
    # Execute steps
    results = {}
    
    # Step 1: Index GWAS files
    if 1 in steps_to_run:
        index_file = step1_index_gwas_files(config, base_dir, args.verbose, args.dry_run)
        results['step1'] = {'index_file': index_file}
    
    # Step 2: Extract GWAS tables
    if 2 in steps_to_run:
        exit_code = step2_extract_gwas_tables(config, base_dir, index_file, args.verbose, args.dry_run)
        results['step2'] = {'exit_code': exit_code}
        if exit_code != 0 and not args.dry_run:
            print("Error in Step 2. Pipeline stopped.")
            sys.exit(exit_code)
    
    # Step 3: Extract significant SNPs
    if 3 in steps_to_run:
        exit_code = step3_extract_significant_snps(config, base_dir, args.verbose, args.dry_run)
        results['step3'] = {'exit_code': exit_code}
        if exit_code != 0 and not args.dry_run:
            print("Error in Step 3. Pipeline stopped.")
            sys.exit(exit_code)
    
    # Step 4: Cross-trait effects
    if 4 in steps_to_run:
        exit_code = step4_cross_trait_effects(config, base_dir, args.verbose, args.dry_run)
        results['step4'] = {'exit_code': exit_code}
        if exit_code != 0 and not args.dry_run:
            print("Error in Step 4. Pipeline stopped.")
            sys.exit(exit_code)
    
    # Step 5: Merge loci
    if 5 in steps_to_run:
        exit_code = step5_merge_loci(config, base_dir, args.verbose, args.dry_run)
        results['step5'] = {'exit_code': exit_code}
        if exit_code != 0 and not args.dry_run:
            print("Error in Step 5. Pipeline stopped.")
            sys.exit(exit_code)
    
    # Step 6: Colocalization
    if 6 in steps_to_run:
        exit_code = step6_run_coloc(config, base_dir, args.verbose, args.dry_run)
        results['step6'] = {'exit_code': exit_code}
        if exit_code != 0 and not args.dry_run:
            print("Error in Step 6. Pipeline stopped.")
            sys.exit(exit_code)
    
    # Step 7: Generate summary
    if 7 in steps_to_run:
        exit_code = step7_generate_summary(config, base_dir, args.verbose, args.dry_run)
        results['step7'] = {'exit_code': exit_code}
        if exit_code != 0 and not args.dry_run:
            print("Error in Step 7. Pipeline stopped.")
            sys.exit(exit_code)
    
    # Final summary
    if args.verbose and not args.dry_run:
        print(f"\n{'='*60}")
        print("Pipeline Execution Summary")
        print("="*60)
        
        successful_steps = [step for step, result in results.items() if 'exit_code' in result and result['exit_code'] == 0]
        
        print(f"Total steps executed: {len(results)}")
        print(f"Successful steps: {len(successful_steps)}")
        
        if len(successful_steps) == len(results):
            print("Status: COMPLETED SUCCESSFULLY")
        else:
            print("Status: COMPLETED WITH ERRORS")
        
        print(f"\nOutput directory: {base_dir}")
        print(f"Results in: {os.path.join(base_dir, config['output']['directories']['snp_analyze'])}")
        
        summary_file = os.path.join(
            base_dir,
            config['output']['directories']['snp_analyze'],
            config['output']['files']['analysis_summary']
        )
        
        if os.path.exists(summary_file):
            print(f"Summary report: {summary_file}")
        
        print(f"\nPipeline finished at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print("="*60)

if __name__ == "__main__":
    main()