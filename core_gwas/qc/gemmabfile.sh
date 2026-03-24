# -k: 指向刚生成的 kinship 文件 (通常在 output/kinship.cXX.txt)
# -c: 指向你的协变量文件 (gwas_ready.cov)
# -p: 指向表型文件 (gwas_ready.pheno)
# -n 1: 跑第一个性状
# -lmm 4: 使用 Wald test (最常用)
# 确保在 workflow 目录下执行
gemma -bfile 01_plink/01_convert/soybean_final \
      -k output/kinship.cXX.txt \
      -c gwas_ready.cov \
      -p gwas_ready.pheno \
      -n 1 -lmm 4 -o gwas_protein_results