# Generalized from soybean project-specific path layout.
#!/usr/bin/env python3
import os
import pandas as pd
import numpy as np

WORKDIR = "data/input/workflow"
os.chdir(WORKDIR)

PHENO_LONG = "../data/phenotype_long.tsv"          # 原始长表
FAM_FILE   = "01_plink/05_final/soybean_wild_final_plink1.fam"  # 103个GDW样本
OUT_DIR    = "03_gwas_input_wild"

TRAITS = ["100SW", "Protein", "Oil"]   # 原始表型列名

print("=== 从 phenotype_long.tsv 重新提取 Wild 大豆表型 ===")

# 1. 读取 wild fam，拿到 103 个 IID（GDWxxx）
fam = pd.read_csv(FAM_FILE, sep=r"\s+", header=None)
fam.columns = ["FID","IID","FA","MO","SEX","PHENO"]
fam["IID"] = fam["IID"].astype(str).str.strip().str.upper()
wild_ids = fam["IID"].unique()
print(f"Wild fam 样本数: {len(wild_ids)}")
print("Wild fam 前5个 IID:", wild_ids[:5])

# 2. 读取原始 phenotype_long.tsv
pheno = pd.read_csv(PHENO_LONG, sep="\t")
pheno.columns = pheno.columns.str.strip()
if "ID" not in pheno.columns:
    raise ValueError("phenotype_long.tsv 中找不到 'ID' 列，请检查列名。")

pheno["ID"] = pheno["ID"].astype(str).str.strip().str.upper()

print("phenotype_long.tsv 总行数:", len(pheno))
print("phenotype_long.tsv 独立ID数:", pheno["ID"].nunique())

# 3. 只保留 ID 在 wild_ids 里的行
pheno_wild = pheno[pheno["ID"].isin(wild_ids)].copy()
print("wild 表型行数:", len(pheno_wild))
print("wild 表型独立ID数:", pheno_wild["ID"].nunique())
print("wild 表型ID示例:", pheno_wild["ID"].unique()[:10])

# 4. 检查是否有缺失列
missing = [t for t in TRAITS if t not in pheno_wild.columns]
if missing:
    raise ValueError(f"phenotype_long.tsv 中缺少这些列: {missing}")

# 5. 计算 Oil:Protein 比例和 log 比例（以原始值为基础，未经标准化）
pheno_wild["Oil_Prot_ratio"] = pheno_wild["Oil"] / pheno_wild["Protein"]
pheno_wild["log_Oil_Prot_ratio"] = np.log(pheno_wild["Oil_Prot_ratio"])

# 6. 按 ID 做多年份均值（不做任何标准化）
traits_all = ["100SW","Protein","Oil","Oil_Prot_ratio","log_Oil_Prot_ratio"]
pheno_wild_avg = pheno_wild.groupby("ID")[traits_all].mean(numeric_only=True).reset_index()
print("wild 平均后独立ID数:", pheno_wild_avg["ID"].nunique())

# 7. 保存一份干净的 wild 表型文件，供后面 prepare_gwas_input_wild 使用
os.makedirs(OUT_DIR, exist_ok=True)
out_pheno_wild = os.path.join(OUT_DIR, "phenotype_wild_avg_with_ratio.tsv")
pheno_wild_avg.to_csv(out_pheno_wild, sep="\t", index=False)
print(f"✓ 已生成: {out_pheno_wild}")
print("预览前5行:")
print(pheno_wild_avg.head())
