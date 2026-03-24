# Generalized from soybean project-specific path layout.
#!/usr/bin/env python3
import os
import pandas as pd

# 路径
PHENO_WILD = "data/input/workflow/03_gwas_input_wild/phenotype_wild_avg_with_ratio.tsv"
FAM_FILE   = "data/input/workflow/01_plink/05_final/soybean_wild_final_plink1.fam"

traits = ["100SW","Protein","Oil","Oil_Prot_ratio","log_Oil_Prot_ratio"]

print("=== 清理 wild 表型中“全缺失”的个体，并同步更新 fam ===")

# 1. 读表型
pheno = pd.read_csv(PHENO_WILD, sep="\t")
pheno["ID"] = pheno["ID"].astype(str).str.strip().str.upper()

print("原始 wild 平均表型 ID 数:", pheno["ID"].nunique())

# 2. 标记“所有性状都 NaN”的行
all_na_mask = pheno[traits].isna().all(axis=1)
to_drop = pheno.loc[all_na_mask, "ID"].unique()
print("需要删除的 ID（所有性状都缺失）个数:", len(to_drop))
print("示例:", to_drop[:10])

# 3. 过滤掉这些 ID
pheno_clean = pheno[~pheno["ID"].isin(to_drop)].copy()
print("清理后剩余 ID 数:", pheno_clean["ID"].nunique())

# 覆盖保存清理后的表型文件
pheno_clean.to_csv(PHENO_WILD, sep="\t", index=False)
print("✓ 已覆盖写回清理后的:", PHENO_WILD)

# 4. 同步更新 wild fam（剔除这些 ID 的基因型）
fam_path = FAM_FILE
fam = pd.read_csv(fam_path, sep=r"\s+", header=None)
fam.columns = ["FID","IID","FA","MO","SEX","PHENO"]
fam["IID"] = fam["IID"].astype(str).str.strip().str.upper()

print("原始 wild fam 行数:", len(fam))

fam_clean = fam[~fam["IID"].isin(to_drop)].copy()
print("清理后 wild fam 行数:", len(fam_clean))

fam_clean.to_csv(fam_path, sep=" ", index=False, header=False)
print("✓ 已覆盖写回清理后的 fam:", fam_path)