# Generalized from soybean project-specific path layout.
# 设置变量
VCF="data/input/data/cleaned/snp.filter.gt.vcf.gz"
PREFIX="01_plink/01_convert/soybean"
QC_PREFIX="01_plink/02_snp_qc/soybean_qc"

# 1. 格式转换 (VCF -> PGEN)
mkdir -p 01_plink/01_convert
plink2 --vcf $VCF \
       --allow-extra-chr \
       --make-pgen \
       --out $PREFIX

# 2. SNP 质量控制 (根据文章推荐参数)
# 文章通常使用 MAF > 0.05, Missing rate < 0.1
mkdir -p 01_plink/02_snp_qc
plink2 --pfile $PREFIX \
       --allow-extra-chr \
       --maf 0.05 \
       --geno 0.1 \
       --snps-only \
       --max-alleles 2 \
       --make-pgen \
       --out $QC_PREFIX