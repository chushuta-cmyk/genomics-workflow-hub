# Generalized from soybean project-specific path layout.
#!/usr/bin/env python3
"""
Haplotype effect analysis for selected lead SNP.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import subprocess
import os
import sys
import logging
from scipy import stats

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Paths
base_dir = "results/cultivated/fine_mapping/log_ratio_17_37608516"
plink_prefix = "data/input/workflow/01_plink/05_final/soybean_final_filtered_bin"
pheno_file = "data/input/workflow/03_gwas_input/phenotype_filtered.tsv"
lead_snp = "17:37608516:A:T"
trait = "log_Oil_Prot_ratio"

# Output files
haplotype_groups_out = os.path.join(base_dir, "haplotype_groups.tsv")
haplotype_effect_out = os.path.join(base_dir, "haplotype_effect_test.tsv")
boxplot_out = os.path.join(base_dir, "haplotype_effect_boxplot.png")

def extract_genotypes():
    """Extract genotypes for lead SNP using plink."""
    # Run plink to export genotypes in raw format
    raw_out = os.path.join(base_dir, "lead_snp_raw")
    cmd = [
        'conda', 'run', '-n', 'gwas', 'plink',
        '--bfile', plink_prefix,
        '--allow-extra-chr',
        '--snp', lead_snp,
        '--recode', 'A',
        '--out', raw_out
    ]
    logger.info(f"Extracting genotypes for SNP {lead_snp}")
    try:
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, check=True)
        logger.info("Plink genotype extraction succeeded.")
    except subprocess.CalledProcessError as e:
        logger.error(f"Plink failed: {e.stderr}")
        # Try direct plink
        cmd2 = ['plink', '--bfile', plink_prefix, '--allow-extra-chr', '--snp', lead_snp, '--recode', 'A', '--out', raw_out]
        try:
            result = subprocess.run(cmd2, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, check=True)
        except subprocess.CalledProcessError as e2:
            logger.error(f"Direct plink also failed: {e2.stderr}")
            sys.exit(1)
    # Load raw file
    raw_file = raw_out + '.raw'
    if not os.path.exists(raw_file):
        logger.error(f"Plink raw file not generated: {raw_file}")
        sys.exit(1)
    df = pd.read_csv(raw_file, sep='\s+')
    # Keep columns: FID, IID, genotype (named like '17:37608516:A:T')
    # The column name may have ':' replaced with '_'? Let's check.
    # Find column that contains the SNP ID
    snp_col = None
    for col in df.columns:
        if lead_snp in col:
            snp_col = col
            break
    if snp_col is None:
        logger.error(f"SNP column not found in raw file.")
        sys.exit(1)
    # Rename for clarity
    df = df[['FID', 'IID', snp_col]].copy()
    df.rename(columns={snp_col: 'dosage'}, inplace=True)
    # dosage is additive coding: 0 = homozygous reference, 1 = heterozygous, 2 = homozygous alternate
    logger.info(f"Genotype extraction complete. Shape: {df.shape}")
    return df

def load_phenotypes():
    """Load phenotype data."""
    df = pd.read_csv(pheno_file, sep='\t')
    # The ID column is 'ID' (GDC***). Need to map to plink IID? 
    # In the .fam file, the individual IDs are also GDC***? Let's check .fam first line.
    # We'll assume IID matches ID column.
    # Keep only ID and trait column
    if trait not in df.columns:
        logger.error(f"Trait column {trait} not found in phenotype file.")
        sys.exit(1)
    # There may be multiple rows per accession (due to years). We need to average across years? 
    # For simplicity, we'll take the mean per accession.
    # Group by ID and compute mean trait value
    pheno_df = df.groupby('ID')[trait].mean().reset_index()
    pheno_df.rename(columns={'ID': 'IID'}, inplace=True)
    logger.info(f"Loaded phenotypes for {len(pheno_df)} individuals.")
    return pheno_df

def merge_and_group(geno_df, pheno_df):
    """Merge genotype and phenotype, assign haplotype groups."""
    merged = pd.merge(geno_df, pheno_df, on='IID', how='inner')
    logger.info(f"Merged data: {len(merged)} samples.")
    # Define groups based on dosage
    # We'll create three groups: 0, 1, 2
    merged['group'] = merged['dosage'].astype(int)
    # Also create allele count grouping: reference allele count (2 - dosage?) Not needed.
    # Save haplotype groups table
    merged.to_csv(haplotype_groups_out, sep='\t', index=False)
    logger.info(f"Saved haplotype groups to {haplotype_groups_out}")
    return merged

def compute_statistics(merged):
    """Compute summary statistics and ANOVA/t-test."""
    # Summary per group
    summary = merged.groupby('group').agg(
        sample_count=('IID', 'count'),
        mean_phenotype=(trait, 'mean'),
        median_phenotype=(trait, 'median'),
        std_phenotype=(trait, 'std'),
        min_phenotype=(trait, 'min'),
        max_phenotype=(trait, 'max')
    ).reset_index()
    # Perform ANOVA across three groups (if at least two groups have >1 sample)
    groups = merged['group'].unique()
    if len(groups) >= 2:
        # One-way ANOVA using scipy
        group_data = [merged.loc[merged['group'] == g, trait] for g in groups]
        f_stat, p_value = stats.f_oneway(*group_data)
        summary['ANOVA_F'] = f_stat
        summary['ANOVA_p'] = p_value
        # Also perform t-test between homozygous groups (0 vs 2) if both exist
        if 0 in groups and 2 in groups:
            t_stat, p_ttest = stats.ttest_ind(
                merged.loc[merged['group'] == 0, trait],
                merged.loc[merged['group'] == 2, trait],
                equal_var=False
            )
            summary['t_test_0vs2_t'] = t_stat
            summary['t_test_0vs2_p'] = p_ttest
        else:
            summary['t_test_0vs2_t'] = np.nan
            summary['t_test_0vs2_p'] = np.nan
    else:
        summary['ANOVA_F'] = np.nan
        summary['ANOVA_p'] = np.nan
        summary['t_test_0vs2_t'] = np.nan
        summary['t_test_0vs2_p'] = np.nan
    # Save effect test table
    summary.to_csv(haplotype_effect_out, sep='\t', index=False)
    logger.info(f"Saved haplotype effect test results to {haplotype_effect_out}")
    return summary

def plot_boxplot(merged):
    """Create boxplot of phenotype by genotype group."""
    fig, ax = plt.subplots(figsize=(8, 6))
    # Boxplot
    box_data = [merged.loc[merged['group'] == g, trait] for g in sorted(merged['group'].unique())]
    ax.boxplot(box_data, labels=sorted(merged['group'].unique()))
    ax.set_xlabel('Genotype group (0=hom ref, 1=het, 2=hom alt)')
    ax.set_ylabel(trait)
    ax.set_title(f'Phenotype distribution by genotype at {lead_snp}')
    # Add sample counts
    for i, g in enumerate(sorted(merged['group'].unique())):
        count = len(merged[merged['group'] == g])
        ax.text(i+1, ax.get_ylim()[0], f'n={count}', ha='center', va='top', fontsize=9)
    plt.tight_layout()
    fig.savefig(boxplot_out, dpi=300)
    logger.info(f"Saved boxplot to {boxplot_out}")
    return fig

def main():
    logger.info("Starting haplotype effect analysis.")
    geno_df = extract_genotypes()
    pheno_df = load_phenotypes()
    merged = merge_and_group(geno_df, pheno_df)
    summary = compute_statistics(merged)
    plot_boxplot(merged)
    logger.info("Haplotype effect analysis completed.")

if __name__ == "__main__":
    main()