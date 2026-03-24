# Generalized from soybean project-specific path layout.
#!/usr/bin/env python3
import os
import pandas as pd

WORKDIR = "data/input/workflow"
os.chdir(WORKDIR)

CLUMP_DIR = "07_clumped_wild"
P_TAG = "1e-064"  # 和上面 out 名里的 P 值标签对应

traits = ["100SW", "Protein", "Oil", "log_ratio"]

print("=== Wild 每个性状 clumping 后的工具数 ===")
for trait in traits:
    clump_path = os.path.join(CLUMP_DIR, f"{trait}.P{P_TAG}.clumped")
    if not os.path.exists(clump_path):
        print(f"{trait}: 找不到 {clump_path}（可能没有显著 SNP）")
        continue
    # PLINK .clumped 是空格分隔、前 2 行 header，字段 SNP 列名通常是 SNP
    df = pd.read_csv(clump_path, delim_whitespace=True, comment="#")
    n_inst = df["SNP"].nunique()
    print(f"{trait}: {n_inst} 个独立工具 SNP")
