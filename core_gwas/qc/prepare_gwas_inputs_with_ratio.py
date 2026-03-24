# Generalized from soybean project-specific path layout.
#!/usr/bin/env python3
import os
import subprocess
import pandas as pd
import numpy as np

WORKDIR = "data/input/workflow"
os.chdir(WORKDIR)

PHENO_FILE   = "../data/phenotype_with_ratio.tsv"
SAMPLES_FILE = "../samples/samples_after_qc.txt"
PFILE_PREFIX = "01_plink/05_final/soybean_final"

TRAITS = ["100SW", "Protein", "Oil", "log_Oil_Prot_ratio"]
ID_COL = "ID"


# 0. 创建目录
os.makedirs("03_gwas_input", exist_ok=True)

# 1. 获取**基因型真实样本列表**（最权威）
print("\n1. 从基因型获取真实样本ID:")
subprocess.run(["plink2", "--pfile", PFILE_PREFIX, "--make-just-fam", "--out", "03_gwas_input/geno_samples"], 
               check=True, stdout=subprocess.DEVNULL)
geno_samples = pd.read_csv("03_gwas_input/geno_samples.fam", sep=r"\s+", header=None)
geno_samples.columns = ["FID", "IID", "FA", "MO", "SEX", "PHENO"]
geno_ids = geno_samples["IID"].astype(str).str.strip().str.upper()
geno_count = len(geno_ids)
print(f"   ✅ 基因型真实样本: {geno_count} 个ID")

# 2. 处理表型
print("\n2. 处理phenotype_with_ratio.tsv:")
pheno = pd.read_csv(PHENO_FILE, sep="\t")
pheno[ID_COL] = pheno[ID_COL].astype(str).str.strip().str.upper()

# 只保留基因型存在的样本，并计算多年份均值
pheno_filt = pheno[pheno[ID_COL].isin(geno_ids)].copy()
pheno_avg = pheno_filt.groupby(ID_COL)[TRAITS].mean(numeric_only=True).reset_index()
pheno_match_count = len(pheno_avg)
print(f"   表型匹配基因型: {pheno_match_count}/{geno_count} 个ID")

# 3. 生成PCA协变量（从eigenvec重建，确保完美对齐）
print("\n3. 重建PCA协变量:")
if os.path.exists("02_pca/pca_data.eigenvec"):
    eigenvec = pd.read_csv("02_pca/pca_data.eigenvec", sep=r"\s+", header=None)
    eigenvec.columns = ["FID", "IID"] + [f"PC{i}" for i in range(1, len(eigenvec.columns)-1)]
else:
    print("   ⚠️ 02_pca/pca_data.eigenvec不存在，跳过PCA，使用仅截距模型")
    eigenvec = pd.DataFrame({"IID": geno_ids, "PC1": 0, "PC2": 0, "PC3": 0})

eigenvec["IID"] = eigenvec["IID"].astype(str).str.strip().str.upper()
pca_clean = eigenvec[eigenvec["IID"].isin(geno_ids)]
print(f"   PCA样本对齐后: {len(pca_clean)} 行")

# 4. 最终合并
sample_order = geno_samples[["IID"]].copy()
merged = sample_order.merge(pca_clean[["IID", "PC1", "PC2", "PC3"]], on="IID", how="left")
merged = merged.merge(pheno_avg, left_on="IID", right_on=ID_COL, how="left")

# --- 新增：Z-score 标准化逻辑 ---
# 我们在填充 -9 之前进行标准化，这样不会把缺失值算进均值里
for trait in TRAITS:
    valid_mask = merged[trait].notna()
    mean_val = merged.loc[valid_mask, trait].mean()
    std_val = merged.loc[valid_mask, trait].std()
    print(f"   标准化 {trait}: mean={mean_val:.3f}, std={std_val:.3f}")
    merged.loc[valid_mask, trait] = (merged.loc[valid_mask, trait] - mean_val) / std_val

# 5. 保存 GEMMA 标准文件
# 先保存一份含 NA 的预览
merged[TRAITS].to_csv("03_gwas_input/gwas_ready.pheno", sep="\t", index=False, header=False, na_rep="NA")

# 处理 GEMMA 专用的 -9 填充
pheno_gemma = merged[TRAITS].fillna(-9) 
pheno_gemma.to_csv("03_gwas_input/gwas_final_clean.pheno", sep="\t", index=False, header=False)

# 协变量保存
merged["Intercept"] = 1
merged[["Intercept", "PC1", "PC2", "PC3"]].to_csv("03_gwas_input/gwas_ready.cov", sep="\t", index=False, header=False)