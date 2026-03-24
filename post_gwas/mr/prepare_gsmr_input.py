# Generalized from soybean project-specific path layout.
#!/usr/bin/env python3
import os
import pandas as pd

WORKDIR = "data/input/workflow"
os.chdir(WORKDIR)

QC_DIR = "05_qc_wild"
OUT_DIR = "09_gsmr_input_wild"
os.makedirs(OUT_DIR, exist_ok=True)

TRAITS = ["100SW", "Protein", "Oil", "log_ratio"]

# wild 群体样本量
N_WILD = 103

print("=== 生成 wild 组 GSMR 输入文件 (SNP A1 A2 freq b se p N) ===")

for trait in TRAITS:
    qc_path = os.path.join(QC_DIR, f"{trait}.qc.tsv")
    if not os.path.exists(qc_path):
        print(f"{trait}: 找不到 {qc_path}，跳过。")
        continue

    df = pd.read_csv(qc_path, sep="\t")

    # 确保必要列存在
    needed = ["SNP", "A1", "A2", "MAF", "BETA", "SE", "P"]
    missing = [c for c in needed if c not in df.columns]
    if missing:
        raise ValueError(f"{qc_path} 缺少列: {missing}")

    out_df = pd.DataFrame({
        "SNP": df["SNP"].astype(str),
        "A1": df["A1"],
        "A2": df["A2"],
        # 假设 GEMMA 的 af 是 A1 的频率（你之前 pipeline 就按这个逻辑用 MAF 了）
        "freq": df["MAF"],
        "b": df["BETA"],
        "se": df["SE"],
        "p": df["P"],
        "N": N_WILD
    })

    out_path = os.path.join(OUT_DIR, f"{trait}.gsmr.txt")
    out_df.to_csv(out_path, sep="\t", index=False)
    print(f"{trait}: GSMR 输入已写出 -> {out_path}")
