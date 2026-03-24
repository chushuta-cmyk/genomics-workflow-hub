# Generalized from soybean project-specific path layout.
#!/usr/bin/env bash
set -euo pipefail

############################
# 1. 激活 conda 环境
############################
source data/miniforge3/etc/profile.d/conda.sh
conda activate gwas

echo "[INFO] Conda env activated:"
which plink2
plink2 --version

############################
# 2. 设置工作目录
############################
WORKDIR=data/input/workflow
cd ${WORKDIR}

echo "[INFO] Working directory: $(pwd)"

############################
# 3. 运行 Snakemake（plink + PCA）
############################
snakemake \
  -s Snakefile.part1_plink_pca \
  --cores 20 \
  --printshellcmds \
  --reason

echo "[INFO] Pipeline finished successfully."
