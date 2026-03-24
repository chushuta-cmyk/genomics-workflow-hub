# Generalized from soybean project-specific path layout.
#!/usr/bin/env bash
set -euo pipefail

cd data/input/workflow

BFILE="01_plink/05_final/soybean_wild_final_plink1"
QC_DIR="05_qc_wild"
LIST_DIR="06_instruments_wild"
CLUMP_DIR="07_clumped_wild"

mkdir -p "${CLUMP_DIR}"

# P 阈值要和 Python 脚本里一致
P_THRESH=1e-6

# trait 名称必须和 qc 文件 & snplist 一致
for TRAIT in 100SW Protein Oil log_ratio; do
  QC_FILE="${QC_DIR}/${TRAIT}.qc.tsv"
  SNPLIST="${LIST_DIR}/${TRAIT}.P${P_THRESH}.snplist"

  if [ ! -f "${SNPLIST}" ]; then
    echo "跳过 ${TRAIT}（没有达到阈值的 SNP 列表）"
    continue
  fi

  echo "=== Clumping ${TRAIT} ==="

  plink \
    --bfile "${BFILE}" \
    --clump "${QC_FILE}" \
    --clump-snp-list "${SNPLIST}" \
    --clump-field P \
    --clump-p1 ${P_THRESH} \
    --clump-p2 1.0 \
    --clump-r2 0.1 \
    --clump-kb 250 \
    --out "${CLUMP_DIR}/${TRAIT}.P${P_THRESH}"

done
