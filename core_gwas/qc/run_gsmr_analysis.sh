# Generalized from soybean project-specific path layout.
#!/bin/bash

################################################################################
# GSMR 因果推断分析 - 完整流程
# 功能: Size → Protein/Oil 的因果关系检验
# 前提: 已完成 convert.sh 格式转换
################################################################################

set -e

# ============================================================================
# 配置参数（根据你的实际路径修改）
# ============================================================================

# 工作目录
WORK_DIR="data/input/workflow"
cd "$WORK_DIR"

# 输入文件
GENO_FILE="01_plink/05_final/soybean_final_plink1"  # PLINK格式基因型
SIZE_GWAS="output/size.ma"                           # Size GWAS汇总数据
PROTEIN_GWAS="output/protein.ma"                     # Protein GWAS汇总数据
OIL_GWAS="output/oil.ma"                             # Oil GWAS汇总数据

# 输出目录
OUT_DIR="05_gsmr_correlation/gsmr_results"
mkdir -p "$OUT_DIR"

# 颜色输出
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m'

echo -e "${BLUE}"
echo "╔════════════════════════════════════════════════════════════════╗"
echo "║       GSMR 因果推断分析 - Size → Protein/Oil                   ║"
echo "╚════════════════════════════════════════════════════════════════╝"
echo -e "${NC}"

# ============================================================================
# 步骤1: 文件检查
# ============================================================================

echo -e "${YELLOW}[步骤1/4] 检查输入文件...${NC}"

# 检查基因型文件
if [[ ! -f "${GENO_FILE}.bed" ]] || [[ ! -f "${GENO_FILE}.bim" ]] || [[ ! -f "${GENO_FILE}.fam" ]]; then
    echo -e "${RED}❌ 错误: 找不到基因型文件 ${GENO_FILE}.{bed,bim,fam}${NC}"
    echo "   当前目录: $(pwd)"
    echo "   请检查路径是否正确"
    exit 1
fi
echo "  ✓ 基因型文件: ${GENO_FILE}"

# 检查GWAS汇总数据
for gwas_file in "$SIZE_GWAS" "$PROTEIN_GWAS" "$OIL_GWAS"; do
    if [[ ! -f "$gwas_file" ]]; then
        echo -e "${RED}❌ 错误: 找不到GWAS文件 $gwas_file${NC}"
        exit 1
    fi
    trait=$(basename "$gwas_file" .ma)
    rows=$(awk 'NR>1' "$gwas_file" | wc -l)
    sig=$(awk 'NR>1 && $7 < 5e-8' "$gwas_file" | wc -l)
    echo "  ✓ $trait: $rows SNPs (p<5e-8: $sig)"
done

# ============================================================================
# 步骤2: GSMR分析 - Size → Protein
# ============================================================================

echo -e "\n${YELLOW}[步骤2/4] GSMR分析: Size → Protein${NC}"

gcta64 --bfile "$GENO_FILE" \
    --gsmr-file <(printf "size\t${SIZE_GWAS}\nprotein\t${PROTEIN_GWAS}\n") \
    --gsmr-direction 2 \
    --effect-plot \
    --gwas-thresh 5e-8 \
    --clump-r2 0.05 \
    --diff-freq 0.2 \
    --gsmr-snp-min 10 \
    --heidi-thresh 0.01 \
    --out "${OUT_DIR}/size_to_protein" 2>&1 | tee "${OUT_DIR}/size_to_protein.log"

if [[ $? -eq 0 ]]; then
    echo -e "${GREEN}✓ Size → Protein 分析完成${NC}"
else
    echo -e "${RED}❌ Size → Protein 分析失败，请检查日志${NC}"
fi

# ============================================================================
# 步骤3: GSMR分析 - Size → Oil
# ============================================================================

echo -e "\n${YELLOW}[步骤3/4] GSMR分析: Size → Oil${NC}"

gcta64 --bfile "$GENO_FILE" \
    --gsmr-file <(printf "size\t${SIZE_GWAS}\noil\t${OIL_GWAS}\n") \
    --gsmr-direction 2 \
    --effect-plot \
    --gwas-thresh 5e-8 \
    --clump-r2 0.05 \
    --diff-freq 0.2 \
    --gsmr-snp-min 10 \
    --heidi-thresh 0.01 \
    --out "${OUT_DIR}/size_to_oil" 2>&1 | tee "${OUT_DIR}/size_to_oil.log"

if [[ $? -eq 0 ]]; then
    echo -e "${GREEN}✓ Size → Oil 分析完成${NC}"
else
    echo -e "${RED}❌ Size → Oil 分析失败，请检查日志${NC}"
fi

# ============================================================================
# 步骤4: 结果汇总
# ============================================================================

echo -e "\n${YELLOW}[步骤4/4] 生成结果报告...${NC}"

python3 << 'PYEOF'
import os
import pandas as pd
from datetime import datetime

out_dir = "05_gsmr_correlation/gsmr_results"

# 读取结果文件
results = []

for analysis in ['size_to_protein', 'size_to_oil']:
    gsmr_file = f"{out_dir}/{analysis}.gsmr"
    
    if os.path.exists(gsmr_file):
        try:
            df = pd.read_csv(gsmr_file, sep='\s+')
            if not df.empty:
                results.append({
                    'Analysis': analysis.replace('_', ' → '),
                    'N_SNPs': df.iloc[0].get('n_SNPs', 'N/A'),
                    'Beta': df.iloc[0].get('bxy', 'N/A'),
                    'SE': df.iloc[0].get('se', 'N/A'),
                    'P_value': df.iloc[0].get('p', 'N/A')
                })
        except Exception as e:
            print(f"⚠️  解析 {analysis} 失败: {e}")
    else:
        print(f"⚠️  未找到结果文件: {gsmr_file}")

# 生成报告
report = []
report.append("=" * 70)
report.append("GSMR 因果推断分析报告")
report.append("=" * 70)
report.append(f"生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

if results:
    report.append("【GSMR结果】")
    report.append("-" * 70)
    for r in results:
        report.append(f"\n{r['Analysis']}:")
        report.append(f"  工具变量SNP数: {r['N_SNPs']}")
        report.append(f"  因果效应 (β):   {r['Beta']}")
        report.append(f"  标准误 (SE):    {r['SE']}")
        report.append(f"  P值:            {r['P_value']}")
        
        try:
            p = float(r['P_value'])
            if p < 0.05:
                report.append(f"  ✓ 显著因果关系 (p < 0.05)")
            else:
                report.append(f"  ⚠️  无显著因果关系 (p ≥ 0.05)")
        except:
            pass
else:
    report.append("⚠️  未找到有效的GSMR结果")

report.append("\n" + "=" * 70)
report.append("\n【输出文件】")
report.append(f"  结果目录: {out_dir}/")
report.append(f"  日志文件: {out_dir}/size_to_*.log")
report.append(f"  GSMR结果: {out_dir}/size_to_*.gsmr")

report_text = "\n".join(report)
print(report_text)

# 保存报告
with open(f"{out_dir}/GSMR_Report.txt", 'w') as f:
    f.write(report_text)

print(f"\n✓ 报告已保存: {out_dir}/GSMR_Report.txt")

PYEOF

# ============================================================================
# 完成
# ============================================================================

echo -e "\n${GREEN}"
echo "╔════════════════════════════════════════════════════════════════╗"
echo "║              ✓ GSMR分析完成!                                   ║"
echo "╚════════════════════════════════════════════════════════════════╝"
echo -e "${NC}"

echo -e "${BLUE}📁 查看结果:${NC}"
echo "   完整报告: cat ${OUT_DIR}/GSMR_Report.txt"
echo "   详细日志: cat ${OUT_DIR}/size_to_protein.log"
echo ""

echo -e "${BLUE}📊 关键文件:${NC}"
ls -lh "${OUT_DIR}"/size_to_*.{gsmr,eff_plot.gz} 2>/dev/null | sed 's/^/   /'
echo ""
