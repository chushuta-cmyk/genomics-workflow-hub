# 运行 Oil (通常是表型文件第 2 列)
gemma -bfile 01_plink/01_convert/soybean_final \
      -k output/kinship.cXX.txt \
      -c gwas_ready.cov \
      -p gwas_ready.pheno \
      -n 2 -lmm 4 -o gwas_oil_results

# 运行 Size (通常是表型文件第 3 列)
gemma -bfile 01_plink/01_convert/soybean_final \
      -k output/kinship.cXX.txt \
      -c gwas_ready.cov \
      -p gwas_ready.pheno \
      -n 3 -lmm 4 -o gwas_size_results