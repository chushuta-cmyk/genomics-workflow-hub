# Generalized from soybean project-specific path layout.
#!/usr/bin/env python3
import os
import pandas as pd

WORKDIR = "data/input/workflow"
os.chdir(WORKDIR)

QC_DIR   = "05_qc_wild"
LIST_DIR = "06_instruments_wild"
os.makedirs(LIST_DIR, exist_ok=True)

TRAITS = ["100SW", "Protein", "Oil", "log_ratio"]

# 设定 P 阈值
P_THRESH = 1e-4 

for trait in TRAITS:
    in_path = os.path.join(QC_DIR, f"{trait}.qc.tsv")
    df = pd.read_csv(in_path, sep="\t")

    # 取显著 SNP
    sig = df[df["P"] < P_THRESH].copy()
    n_sig = len(sig)
    print(f"{trait}: P<{P_THRESH:g} 候选 SNP 数 = {n_sig}")

    if n_sig == 0:
        print(f"  -> 暂无显著 SNP，后面 clumping 会跳过。")
        continue

    # 为 PLINK 准备一个简单列表：SNP ID 一列
    list_path = os.path.join(LIST_DIR, f"{trait}.P{P_THRESH:g}.snplist")
    sig["SNP"].to_csv(list_path, index=False, header=False)
    print("  -> SNP 列表已保存:", list_path)
