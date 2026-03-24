#!/usr/bin/env python3
"""
Fine mapping for selected GWAS peak.
Input: cultivated GWAS summary for log_ratio, genotype data, gene annotation.
Output: regional plots, LD heatmap, candidate genes.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import subprocess
import os
import sys
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Paths
base_dir = "/data03/karama/projects/soybean_analysis/cultivated/fine_mapping/log_ratio_17_37608516"
gwas_file = "/data03/karama/projects/soybean_analysis/cultivated/gwas/gwas_log_ratio.tsv"
plink_prefix = "/data03/karama/soybean_gwas_filtered/workflow/01_plink/05_final/soybean_final_filtered_bin"
gtf_file = "/data03/karama/soybean_gwas/reference/W82.a6.v1/annotation/Gmax_880_Wm82.a6.v1.gene_exons.gtf"
locus_summary_file = "/data03/karama/projects/soybean_analysis/cultivated/post_gwas/locus_summary_log_ratio.tsv"

# Selected locus info
selected_chr = "17"
selected_pos = 37608516
selected_lead_snp = f"{selected_chr}:{selected_pos}:A:T"  # from sig_all_traits
window = 250000  # base pairs
# Chromosome naming in GTF
gtf_chr = f"Gm{selected_chr}"

# Create output directory
os.makedirs(base_dir, exist_ok=True)

def load_locus_metadata():
    """Load locus metadata from locus summary file."""
    df = pd.read_csv(locus_summary_file, sep='\t')
    # Find row matching lead SNP position
    mask = df['lead_snp_pos'] == selected_pos
    if not mask.any():
        # try matching lead_snp string
        mask = df['lead_snp'] == selected_lead_snp
    if not mask.any():
        logger.error(f"Locus with lead SNP {selected_lead_snp} not found in {locus_summary_file}")
        sys.exit(1)
    locus = df[mask].iloc[0]
    logger.info(f"Found locus: {locus['locus_id']}")
    return locus

def create_peak_metadata(locus):
    """Create peak_metadata.tsv."""
    metadata = {
        'trait': 'log_ratio',
        'locus_id': locus['locus_id'],
        'chr': locus['chr'],
        'start': locus['start'],
        'end': locus['end'],
        'locus_size': locus['locus_size'],
        'num_snps': locus['num_snps'],
        'lead_snp': locus['lead_snp'],
        'lead_snp_pos': locus['lead_snp_pos'],
        'lead_p': locus['lead_p'],
        'traits_in_locus': locus['traits_in_locus'],
        'locus_width': locus['locus_width'],
        'fine_mapping_window_start': selected_pos - window,
        'fine_mapping_window_end': selected_pos + window,
        'fine_mapping_window_bp': window * 2
    }
    df = pd.DataFrame([metadata])
    out_path = os.path.join(base_dir, "peak_metadata.tsv")
    df.to_csv(out_path, sep='\t', index=False)
    logger.info(f"Saved peak metadata to {out_path}")
    return df

def extract_local_gwas():
    """Extract GWAS SNPs within window."""
    logger.info(f"Loading GWAS file {gwas_file}")
    # Determine column names
    with open(gwas_file, 'r') as f:
        header = f.readline().strip().split('\t')
    cols = ['chr', 'rs', 'ps', 'allele1', 'allele0', 'af', 'beta', 'se', 'p_wald']
    missing = [c for c in cols if c not in header]
    if missing:
        logger.warning(f"Columns {missing} not found. Using available columns.")
        usecols = None
    else:
        usecols = cols
    df = pd.read_csv(gwas_file, sep='\t', usecols=usecols, low_memory=False)
    # Filter chromosome and position
    df['chr'] = df['chr'].astype(str)
    region_df = df[(df['chr'] == selected_chr) & 
                   (df['ps'] >= selected_pos - window) & 
                   (df['ps'] <= selected_pos + window)].copy()
    # Sort by position
    region_df = region_df.sort_values('ps')
    out_path = os.path.join(base_dir, "local_gwas_region.tsv")
    region_df.to_csv(out_path, sep='\t', index=False)
    logger.info(f"Extracted {len(region_df)} SNPs in region. Saved to {out_path}")
    return region_df

def plot_local_gwas(region_df):
    """Create local GWAS zoom Manhattan plot."""
    fig, ax = plt.subplots(figsize=(10, 6))
    region_df['log10p'] = -np.log10(region_df['p_wald'])
    lead_mask = region_df['ps'] == selected_pos
    ax.scatter(region_df['ps'], region_df['log10p'], s=20, alpha=0.6, c='steelblue', edgecolors='none')
    if lead_mask.any():
        ax.scatter(region_df.loc[lead_mask, 'ps'], region_df.loc[lead_mask, 'log10p'], 
                  s=100, c='red', edgecolors='black', label='Lead SNP')
    ax.axhline(y=-np.log10(5e-8), color='r', linestyle='--', alpha=0.5, label='p=5e-8')
    ax.set_xlabel(f'Chromosome {selected_chr} position (bp)')
    ax.set_ylabel('-log10(p)')
    ax.set_title(f'Regional association plot for log_ratio\nchr{selected_chr}:{selected_pos - window}-{selected_pos + window}')
    ax.legend()
    plt.tight_layout()
    plot_path = os.path.join(base_dir, "local_gwas_zoom.png")
    fig.savefig(plot_path, dpi=300)
    logger.info(f"Saved local GWAS zoom plot to {plot_path}")
    return fig

def compute_ld_matrix(region_df):
    """Compute pairwise LD matrix using plink (via conda gwas environment)."""
    snp_ids = region_df['rs'].tolist()
    snp_list_file = os.path.join(base_dir, "snp_list.txt")
    with open(snp_list_file, 'w') as f:
        for snp in snp_ids:
            f.write(snp + '\n')
    ld_out_prefix = os.path.join(base_dir, "local_ld")
    # Use conda run to invoke plink from gwas environment
    cmd = [
        'conda', 'run', '-n', 'gwas', 'plink',
        '--bfile', plink_prefix,
        '--allow-extra-chr',
        '--chr', selected_chr,
        '--extract', snp_list_file,
        '--r2', 'square',
        '--out', ld_out_prefix
    ]
    logger.info(f"Running plink LD calculation...")
    try:
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, check=True)
        logger.info(f"Plink succeeded.")
    except subprocess.CalledProcessError as e:
        logger.error(f"Plink failed: {e.stderr}")
        # Try direct plink (maybe in PATH)
        cmd2 = ['plink', '--bfile', plink_prefix, '--allow-extra-chr', '--chr', selected_chr,
                '--extract', snp_list_file, '--r2', 'square', '--out', ld_out_prefix]
        try:
            result = subprocess.run(cmd2, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, check=True)
            logger.info("Direct plink succeeded.")
        except subprocess.CalledProcessError as e2:
            logger.error(f"Direct plink also failed: {e2.stderr}")
            return None
    # Load LD matrix
    ld_file = ld_out_prefix + '.ld'
    if os.path.exists(ld_file):
        ld_matrix = pd.read_csv(ld_file, sep='\s+', header=None)
        ld_matrix.to_csv(os.path.join(base_dir, "local_ld_matrix.tsv"), sep='\t', index=False, header=False)
        logger.info(f"Saved LD matrix to local_ld_matrix.tsv")
        return ld_matrix
    else:
        logger.warning(f"Plink did not produce .ld file")
        return None

def plot_ld_heatmap(ld_matrix, region_df):
    """Create LD heatmap from pairwise r²."""
    if ld_matrix is None:
        logger.warning("No LD matrix, skipping heatmap.")
        return
    fig, ax = plt.subplots(figsize=(8, 8))
    sns.heatmap(ld_matrix, cmap='RdBu_r', center=0, square=True, 
                xticklabels=False, yticklabels=False, ax=ax)
    ax.set_title('Pairwise LD (r²) heatmap')
    # Add chromosome position ticks for a subset of SNPs
    step = max(1, len(region_df) // 10)
    tick_positions = np.arange(0, len(region_df), step)
    tick_labels = region_df.iloc[tick_positions]['ps'].astype(str).tolist()
    ax.set_xticks(tick_positions)
    ax.set_yticks(tick_positions)
    ax.set_xticklabels(tick_labels, rotation=90, fontsize=8)
    ax.set_yticklabels(tick_labels, fontsize=8)
    plt.tight_layout()
    heatmap_path = os.path.join(base_dir, "local_ld_heatmap.png")
    fig.savefig(heatmap_path, dpi=300)
    logger.info(f"Saved LD heatmap to {heatmap_path}")
    return fig

def extract_genes():
    """Extract genes within window from GTF."""
    logger.info(f"Loading gene annotation from {gtf_file}")
    genes = {}
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            feature = parts[2]
            # consider both 'gene' and 'transcript' as gene-level features
            if feature not in ['gene', 'transcript']:
                continue
            chrom = parts[0]
            if chrom != gtf_chr:
                continue
            start = int(parts[3])
            end = int(parts[4])
            if start > selected_pos + window or end < selected_pos - window:
                continue
            attrs = parts[8]
            gene_id = None
            gene_name = None
            for attr in attrs.split(';'):
                attr = attr.strip()
                if attr.startswith('gene_id'):
                    gene_id = attr.split('"')[1] if '"' in attr else attr.split(' ')[1]
                elif attr.startswith('gene_name'):
                    gene_name = attr.split('"')[1] if '"' in attr else attr.split(' ')[1]
                elif attr.startswith('Name'):
                    gene_name = attr.split('"')[1] if '"' in attr else attr.split(' ')[1]
            if gene_id is None:
                continue
            # deduplicate by gene_id, keep the one with smallest start
            if gene_id not in genes or start < genes[gene_id]['start']:
                genes[gene_id] = {
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'strand': parts[6],
                    'gene_id': gene_id,
                    'gene_name': gene_name,
                    'attributes': attrs
                }
    genes_df = pd.DataFrame(list(genes.values()))
    out_path = os.path.join(base_dir, "candidate_genes_in_window.tsv")
    genes_df.to_csv(out_path, sep='\t', index=False)
    logger.info(f"Found {len(genes_df)} genes in window. Saved to {out_path}")
    return genes_df

def main():
    logger.info("Starting fine mapping for selected locus.")
    # Step 1: Load locus metadata
    locus = load_locus_metadata()
    # Step 2: Create peak metadata
    create_peak_metadata(locus)
    # Step 3: Extract local GWAS region
    region_df = extract_local_gwas()
    # Step 4: Plot local GWAS zoom
    plot_local_gwas(region_df)
    # Step 5: Compute LD matrix
    ld_matrix = compute_ld_matrix(region_df)
    # Step 6: Plot LD heatmap
    plot_ld_heatmap(ld_matrix, region_df)
    # Step 7: Extract genes
    genes_df = extract_genes()
    logger.info("Fine mapping steps completed.")

if __name__ == "__main__":
    main()