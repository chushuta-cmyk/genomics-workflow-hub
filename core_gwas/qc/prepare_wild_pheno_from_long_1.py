# Generalized from soybean project-specific path layout.
#!/usr/bin/env python3
import os
import pandas as pd
import numpy as np

# 绝对路径
PHENO_LONG = "data/input/data/phenotype_original.tsv"
FAM_FILE   = "data/input/workflow/01_plink/05_final/soybean_wild_final_plink1.fam"
OUT_DIR    = "data/input/workflow/03_gwas_input_wild"

TRAITS = ["100SW", "Protein", "Oil"]

print("=== 从 phenotype_original.tsv 重新提取 Wild 大豆表型 ===")

# 1. 读取 wild fam
fam = pd.read_csv(FAM_FILE, sep=r"\s+", header=None)
fam.columns = ["FID","IID","FA","MO","SEX","PHENO"]
fam["IID"] = fam["IID"].astype(str).str.strip().str.upper()
wild_ids = fam["IID"].unique()
print(f"Wild fam 样本数: {len(wild_ids)}")
print("Wild fam 前5个 IID:", wild_ids[:5])

# 2. 读取原始 phenotype_original.tsv
pheno = pd.read_csv(PHENO_LONG, sep="\t")
pheno.columns = pheno.columns.str.strip()
if "ID" not in pheno.columns:
    raise ValueError("phenotype_original.tsv 中找不到 'ID' 列，请检查列名。")

pheno["ID"] = pheno["ID"].astype(str).str.strip().str.upper()

print("phenotype_original.tsv 总行数:", len(pheno))
print("phenotype_original.tsv 独立ID数:", pheno["ID"].nunique())

# 3. 过滤出 wild ID
pheno_wild = pheno[pheno["ID"].isin(wild_ids)].copy()
print("wild 表型行数:", len(pheno_wild))
print("wild 表型独立ID数:", pheno_wild["ID"].nunique())
print("wild 表型ID示例:", pheno_wild["ID"].unique()[:10])

missing = [t for t in TRAITS if t not in pheno_wild.columns]
if missing:
    raise ValueError(f"phenotype_original.tsv 中缺少这些列: {missing}")

# 4. 计算 Oil:Protein 比例（原始）和 log 比例
pheno_wild["Oil_Prot_ratio"] = pheno_wild["Oil"] / pheno_wild["Protein"]
pheno_wild["log_Oil_Prot_ratio"] = np.log(pheno_wild["Oil_Prot_ratio"])

# 5. 按 ID 做均值（不做标准化）
traits_all = ["100SW","Protein","Oil","Oil_Prot_ratio","log_Oil_Prot_ratio"]
pheno_wild_avg = pheno_wild.groupby("ID")[traits_all].mean(numeric_only=True).reset_index()
print("wild 平均后独立ID数:", pheno_wild_avg["ID"].nunique())

# 6. 保存
os.makedirs(OUT_DIR, exist_ok=True)
out_pheno_wild = os.path.join(OUT_DIR, "phenotype_wild_avg_with_ratio.tsv")
pheno_wild_avg.to_csv(out_pheno_wild, sep="\t", index=False)
print(f"✓ 已生成: {out_pheno_wild}")
print("预览前5行:")
print(pheno_wild_avg.head())
