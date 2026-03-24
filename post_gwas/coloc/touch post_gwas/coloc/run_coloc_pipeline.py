#!/usr/bin/env python3

"""
Colocalization pipeline
Input: cross-trait GWAS results
Output: coloc summary
"""

import os

INPUT = "results/post_gwas/cross_trait.tsv"
OUTPUT = "results/post_gwas/coloc"

os.makedirs(OUTPUT, exist_ok=True)


def run(cmd):
    print(f"\n[RUN] {cmd}\n")
    os.system(cmd)


# =========================
# STEP 1: Convert format
# =========================
run(f"python post_gwas/coloc/ma_convert.py "
    f"--input {INPUT} "
    f"--output {OUTPUT}/coloc_input.tsv")

# =========================
# STEP 2: Run coloc
# =========================
run(f"Rscript post_gwas/coloc/run_coloc.R "
    f"{OUTPUT}/coloc_input.tsv "
    f"{OUTPUT}/coloc_raw.tsv")

# =========================
# STEP 3: Summarize
# =========================
run(f"python post_gwas/coloc/generate_coloc_summary.py "
    f"--input {OUTPUT}/coloc_raw.tsv "
    f"--output {OUTPUT}/coloc_summary.tsv")

print("✅ coloc pipeline done")