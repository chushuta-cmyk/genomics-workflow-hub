#!/usr/bin/env python3

"""
Post-GWAS unified pipeline
DO NOT RUN GWAS again
Only process existing summary statistics
"""

import os

# ===== INPUT =====
GWAS_DIR = "data/gwas_results"   # 通用路径
OUTPUT_DIR = "results/post_gwas"

os.makedirs(OUTPUT_DIR, exist_ok=True)


def run(cmd):
    print(f"\n[RUN] {cmd}\n")
    os.system(cmd)


# ========================
# STEP 1: Extract significant SNPs
# ========================
run(f"python post_gwas/significant_snps/extract_significant_snps.py "
    f"--input {GWAS_DIR} "
    f"--output {OUTPUT_DIR}/significant_snps.tsv")

# ========================
# STEP 2: Merge loci
# ========================
run(f"python post_gwas/locus_merging/merge_loci.py "
    f"--input {OUTPUT_DIR}/significant_snps.tsv "
    f"--window 250000 "
    f"--output {OUTPUT_DIR}/loci.tsv")

# ========================
# STEP 3: Cross-trait analysis
# ========================
run(f"python post_gwas/cross_trait/cross_trait_effects.py "
    f"--input {OUTPUT_DIR}/loci.tsv "
    f"--output {OUTPUT_DIR}/cross_trait.tsv")

# ========================
# STEP 4: Colocalization
# ========================
run(f"Rscript post_gwas/coloc/run_coloc.R "
    f"{OUTPUT_DIR}/cross_trait.tsv "
    f"{OUTPUT_DIR}/coloc_results.tsv")

# ========================
# STEP 5: Comparison
# ========================
run(f"python post_gwas/comparison/compare_gwas.py "
    f"--input {OUTPUT_DIR}/loci.tsv "
    f"--output {OUTPUT_DIR}/comparison.tsv")

print("\n✅ Post-GWAS pipeline completed.")