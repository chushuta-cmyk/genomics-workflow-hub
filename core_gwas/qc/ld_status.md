# LD Analysis Status (Wild Soybean)

## LD Decay Calculation
- **Status:** Completed
- **Command:** `plink --bfile ... --ld-window 99999 --ld-window-kb 1000 --ld-window-r2 0 --r2 --out wild_ld_decay_full`
- **Output file:** `ld_analysis/wild_ld_decay_full.ld` (≈18.7 GB)
- **Log file:** `ld_analysis/wild_ld_decay_full.log`
- **Completion time:** 2026‑03‑09 15:42:12
- **Variants processed:** 428,454
- **Samples:** 103 wild soybean accessions

**Missing downstream outputs:**
- Binned LD decay summary (mean \(r^2\) vs. distance)
- LD decay curve plot
- Average \(r^2\) at given distance thresholds (e.g., \(r^2 < 0.2\))

## Haplotype Block Detection
- **Status:** Not completed (phenotype requirement)
- **Command:** `plink --bfile ... --blocks`
- **Log file:** `ld_analysis/wild_haplotype_blocks.log`
- **Issue:** PLINK requires at least two founders with non‑missing phenotypes. The default can be overridden with `--blocks no-pheno-req`.
- **Raw haplotype block output:** None

## Recommended Actions
1. Run an LD decay summary script (e.g., `scripts/ld_decay_summary.py`) on the `.ld` file.
2. Re‑run haplotype block detection with `--blocks no-pheno-req`.
3. Overlap GWAS loci with LD blocks for fine‑mapping.

---
*Generated on 2026‑03‑09*