# Generalized from soybean project-specific path layout.
#!/usr/bin/env python3
import pandas as pd
import os
from pathlib import Path

project_root = Path("results/wild")
post_gwas_dir = project_root / "post_gwas"
reports_dir = project_root / "reports"

# Read significant SNP summary
sig_summary_path = post_gwas_dir / "significant_snps" / "sig_summary.txt"
with open(sig_summary_path) as f:
    lines = f.readlines()
# parse lines
sig_counts = {}
for line in lines:
    if ':' in line and 'SNPs' in line and line.strip().endswith('SNPs'):
        parts = line.strip().split(':')
        trait = parts[0].strip()
        count = int(parts[1].strip().split()[0])
        sig_counts[trait] = count

# Read locus summary stats for each trait
traits = ['100SW', 'Protein', 'Oil', 'log_ratio']
locus_stats = {}
for trait in traits:
    stats_file = post_gwas_dir / f"locus_summary_{trait}.stats.tsv"
    if stats_file.exists():
        df = pd.read_csv(stats_file, sep='\t')
        # convert to dict
        stats = dict(zip(df['metric'], df['value']))
        locus_stats[trait] = stats
    else:
        locus_stats[trait] = {}

# Read cross-trait effect summary
cross_summary_file = post_gwas_dir / "cross_trait_effect.summary.tsv"
cross_df = pd.read_csv(cross_summary_file, sep='\t')
cross_summary = dict(zip(cross_df['significance_pattern'], cross_df['count']))

# Build final summary table
rows = []
for trait in traits:
    row = {
        'trait': trait,
        'significant_snps': sig_counts.get(trait, 0),
        'total_loci': locus_stats.get(trait, {}).get('total_loci', 0),
        'mean_locus_size': locus_stats.get(trait, {}).get('mean_locus_size', 0),
        'median_locus_size': locus_stats.get(trait, {}).get('median_locus_size', 0),
        'max_locus_size': locus_stats.get(trait, {}).get('max_locus_size', 0),
        'mean_snps_per_locus': locus_stats.get(trait, {}).get('mean_snps_per_locus', 0),
        'single_snp_loci': locus_stats.get(trait, {}).get('single_snp_loci', 0),
        'multi_trait_loci': locus_stats.get(trait, {}).get('multi_trait_loci', 0),
    }
    rows.append(row)

# Overall row (aggregate across traits)
overall_row = {
    'trait': 'OVERALL',
    'significant_snps': sum(sig_counts.values()),
    'total_loci': sum(locus_stats[t].get('total_loci', 0) for t in traits),
    'mean_locus_size': sum(locus_stats[t].get('mean_locus_size', 0) for t in traits) / len(traits),
    'median_locus_size': sum(locus_stats[t].get('median_locus_size', 0) for t in traits) / len(traits),
    'max_locus_size': max(locus_stats[t].get('max_locus_size', 0) for t in traits),
    'mean_snps_per_locus': sum(locus_stats[t].get('mean_snps_per_locus', 0) for t in traits) / len(traits),
    'single_snp_loci': sum(locus_stats[t].get('single_snp_loci', 0) for t in traits),
    'multi_trait_loci': sum(locus_stats[t].get('multi_trait_loci', 0) for t in traits),
}
rows.append(overall_row)

df_summary = pd.DataFrame(rows)
output_path = reports_dir / "final_merged_summary.tsv"
df_summary.to_csv(output_path, sep='\t', index=False)
print(f"Saved final merged summary to {output_path}")

# Also create a simple text summary for the report
with open(reports_dir / "final_summary.txt", 'w') as f:
    f.write("Wild Soybean GWAS Summary\n")
    f.write("==========================\n\n")
    f.write("Significant SNPs per trait (p < 5e-08):\n")
    for trait, count in sig_counts.items():
        f.write(f"  {trait}: {count} SNPs\n")
    f.write(f"  Total unique SNPs: {overall_row['significant_snps']}\n\n")
    f.write("Locus statistics per trait:\n")
    for trait in traits:
        stats = locus_stats.get(trait, {})
        f.write(f"  {trait}:\n")
        f.write(f"    Total loci: {stats.get('total_loci', 0)}\n")
        f.write(f"    Mean locus size: {stats.get('mean_locus_size', 0):.2f} bp\n")
        f.write(f"    Max locus size: {stats.get('max_locus_size', 0)} bp\n")
        f.write(f"    Single SNP loci: {stats.get('single_snp_loci', 0)}\n")
        f.write(f"    Multi‑trait loci: {stats.get('multi_trait_loci', 0)}\n")
    f.write("\nCross‑trait patterns (SNPs significant in multiple traits):\n")
    for pattern, count in cross_summary.items():
        f.write(f"  {pattern}: {count} SNPs\n")

print("Created final_summary.txt")