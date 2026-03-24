# Generalized from soybean project-specific path layout.
#!/bin/bash
set -e
source data/miniforge3/etc/profile.d/conda.sh
conda activate gwas

plink \
--bfile data/input/workflow/01_plink/05_final/soybean_final_filtered_bin \
--allow-extra-chr \
--blocks no-pheno-req \
--out cultivated_haplotype_blocks

echo "Haplotype block detection completed for cultivated soybean."