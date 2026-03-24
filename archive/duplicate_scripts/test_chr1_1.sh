#!/bin/bash
set -e
source /data03/karama/miniforge3/etc/profile.d/conda.sh
conda activate gwas

plink \
--bfile /data03/karama/soybean_gwas_filtered/workflow/01_plink/05_final/soybean_final_filtered_bin \
--chr 1 \
--allow-extra-chr \
--r2 \
--ld-window-kb 1000 \
--ld-window 99999 \
--ld-window-r2 0 \
--out cultivated_ld_decay_chr1

echo "Test LD calculation for chromosome 1 completed for cultivated soybean."