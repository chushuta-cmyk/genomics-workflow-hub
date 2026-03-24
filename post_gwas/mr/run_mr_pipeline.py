#!/usr/bin/env python3

"""
Mendelian Randomization Pipeline
Input: GWAS summary stats
Output: MR results
"""

import os

INPUT = "results/post_gwas/significant_snps.tsv"
OUTPUT = "results/post_gwas/mr"

os.makedirs(OUTPUT, exist_ok=True)


def run(cmd):
    print(f"\n[RUN] {cmd}\n")
    os.system(cmd)


# =========================
# STEP 1: Instrument selection (clumping)
# =========================
run(f"python post_gwas/mr/select_instruments.py "
    f"--input {INPUT} "
    f"--output {OUTPUT}/instruments.tsv")

# =========================
# STEP 2: Harmonization
# =========================
run(f"python post_gwas/mr/harmonize_data.py "
    f"--input {OUTPUT}/instruments.tsv "
    f"--output {OUTPUT}/harmonized.tsv")

# =========================
# STEP 3: Prepare GSMR input
# =========================
run(f"python post_gwas/mr/prepare_gsmr_input.py "
    f"--input {OUTPUT}/harmonized.tsv "
    f"--output {OUTPUT}/gsmr_input.tsv")

# =========================
# STEP 4: Run MR (IVW / Egger)
# =========================
run(f"python post_gwas/mr/run_mr_ivw_egger.py "
    f"--input {OUTPUT}/gsmr_input.tsv "
    f"--output {OUTPUT}/mr_results.tsv")

# =========================
# STEP 5: Plot
# =========================
run(f"python post_gwas/mr/plot_mr_results.py "
    f"--input {OUTPUT}/mr_results.tsv "
    f"--output {OUTPUT}/plots/")

print("✅ MR pipeline completed.")