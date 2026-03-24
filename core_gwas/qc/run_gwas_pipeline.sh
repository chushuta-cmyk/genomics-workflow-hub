#!/bin/bash

# --- 参数配置 ---
PHENOTYPE="phenotype_long.tsv"
# 【请根据你的 ls 结果修改这里】
PCA_VEC="01_plink/plink.eigenvec" 
GENOTYPE_PREFIX="01_plink/your_genotype_prefix" # 你的 .bed 文件前缀

echo ">>> 步骤 1: 预处理表型数据与协变量..."
if [ ! -f "$PCA_VEC" ]; then
    echo "错误: 找不到 PCA 结果文件 $PCA_VEC，请检查路径。"
    exit 1
fi
# 运行 Python (确保 Python 脚本里的路径也同步更新或通过参数传递)
python3 prepare_gwas_input.py

echo ">>> 步骤 2: 运行 GWAS (GEMMA)..."
if ! command -v gemma &> /dev/null; then
    echo "错误: 未找到 gemma，请运行 'conda install -c bioconda gemma' 安装。"
    exit 1
fi

# 计算 Kinship
gemma -bfile $GENOTYPE_PREFIX -gk 1 -o kinship
# 运行 LMM (假设跑第 1 个性状 Protein)
gemma -bfile $GENOTYPE_PREFIX -k output/kinship.cXX.txt \
      -c gwas_ready.cov -p gwas_ready.pheno -n 1 -lmm 4 -o gwas_res

echo ">>> 步骤 3: 绘制曼哈顿图与 Q-Q 图..."
if [ -f "analysis_plots.R" ]; then
    Rscript analysis_plots.R GWAS output/gwas_res.assoc.txt
else
    echo "错误: 找不到 analysis_plots.R 脚本。"
fi