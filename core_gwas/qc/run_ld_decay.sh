# Generalized from soybean project-specific path layout.
#!/bin/bash
set -e
source data/miniforge3/etc/profile.d/conda.sh
conda activate gwas

plink \
--bfile data/input/workflow/01_plink/05_final/soybean_wild_chr1_20 \
--allow-extra-chr \
--r2 \
--ld-window-kb 1000 \
--ld-window 99999 \
--ld-window-r2 0 \
--out wild_ld_decay

echo "LD decay calculation completed."