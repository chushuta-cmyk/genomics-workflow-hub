# Generalized from soybean project-specific path layout.
#!/bin/bash
set -e
source data/miniforge3/etc/profile.d/conda.sh
conda activate gwas

plink \
--bfile data/input/workflow/01_plink/05_final/soybean_wild_chr1_20 \
--allow-extra-chr \
--blocks \
--out wild_haplotype_blocks

echo "Haplotype block detection completed."
