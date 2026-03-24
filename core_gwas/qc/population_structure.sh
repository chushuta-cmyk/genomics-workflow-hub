# 1. LD Pruning with Unique IDs
plink2 --pfile 01_plink/02_snp_qc/soybean_qc \
       --allow-extra-chr \
       --set-all-var-ids '@:#:$r:$a' \
       --indep-pairwise 50 10 0.2 \
       --out 01_plink/03_ld_prune/soybean_pruned

# 2. PCA using the pruned list
plink2 --pfile 01_plink/02_snp_qc/soybean_qc \
       --allow-extra-chr \
       --set-all-var-ids '@:#:$r:$a' \
       --extract 01_plink/03_ld_prune/soybean_pruned.prune.in \
       --pca 10 \
       --out 01_plink/04_pca/soybean_pca