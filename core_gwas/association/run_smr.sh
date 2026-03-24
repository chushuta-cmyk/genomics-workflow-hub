# SMR 软件通常需要指定 GWAS 的 .ma 文件和 eQTL 的前缀
smr --nmr-summary protein_gwas.ma \
    --beql-summary 02_eqtl/soybean_bean_eqtl \
    --out protein_smr_res \
    --thread 8 \
    --peqtl-smr 5e-5 \
    --diff-freq-prop 0.1