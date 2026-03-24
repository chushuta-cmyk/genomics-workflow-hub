# Generalized from soybean project-specific path layout.
#!/usr/bin/env bash
set -euo pipefail

cd data/input/workflow

# wild 的 LD 参考面板
BFILE="01_plink/05_final/soybean_wild_chr1_20"

GSMR_DIR="10_gsmr_wild"
mkdir -p "${GSMR_DIR}"

# GSMR 输入（上一步生成）
IN_DIR="09_gsmr_input_wild"

# 统一参数
P_GSMR=1e-4
HEIDI_THRESH=0.01
SNP_MIN=10

run_gsmr_pair () {
  local EXP=$1
  local OUT=$2

  local EXP_FILE="${IN_DIR}/${EXP}.gsmr.txt"
  local OUT_FILE="${IN_DIR}/${OUT}.gsmr.txt"
  local OUT_PREFIX="${GSMR_DIR}/${EXP}_to_${OUT}"

  echo "=== GSMR: ${EXP} -> ${OUT} ==="
  echo "  暴露文件:  ${EXP_FILE}"
  echo "  结局文件:  ${OUT_FILE}"
  echo "  输出前缀:  ${OUT_PREFIX}"

  gcta64 \
    --bfile "${BFILE}" \
    --gsmr-file "${EXP_FILE}" "${OUT_FILE}" \
    --gsmr-direction 0 \
    --gsmr-p "${P_GSMR}" \
    --gsmr-snp-min "${SNP_MIN}" \
    --heidi-thresh "${HEIDI_THRESH}" \
    --out "${OUT_PREFIX}"
}

# size(100SW) ↔ Oil
run_gsmr_pair 100SW Oil
run_gsmr_pair Oil 100SW

# size ↔ Protein
run_gsmr_pair 100SW Protein
run_gsmr_pair Protein 100SW

# size ↔ log_ratio
run_gsmr_pair 100SW log_ratio
run_gsmr_pair log_ratio 100SW

# Oil ↔ Protein
run_gsmr_pair Oil Protein
run_gsmr_pair Protein Oil
