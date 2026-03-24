# Generalized from soybean project-specific path layout.
#!/usr/bin/env python3
import os
import pandas as pd
import numpy as np

WORKDIR = "data/input/workflow"
os.chdir(WORKDIR)

IN_DIR  = "04_gemma_wild"
OUT_DIR = "05_qc_wild"
os.makedirs(OUT_DIR, exist_ok=True)

# wild 四个主要性状：1=100SW, 2=Protein, 3=Oil, 5=log_ratio
TRAITS = {
    "100SW":  "wild_trait_1.assoc.txt",
    "Protein": "wild_trait_2.assoc.txt",
    "Oil":     "wild_trait_3.assoc.txt",
    "log_ratio": "wild_trait_5.assoc.txt",
}

# MAF 阈值，可以以后根据情况改
MAF_MIN = 0.01

def load_and_standardize(path):
    df = pd.read_csv(path, delim_whitespace=True, low_memory=False)

    # 标准列名
    rename_map = {
        "chr": "CHR",
        "rs": "SNP",
        "ps": "BP",
        "allele1": "A1",
        "allele0": "A2",
        "af": "MAF",
        "beta": "BETA",
        "se": "SE",
        "p_wald": "P"
    }
    df = df.rename(columns=rename_map)

    # 只保留我们关心的列（如果有别的质量指标，后面可以加上来）
    keep_cols = [c for c in ["SNP","CHR","BP","A1","A2","MAF","BETA","SE","P"] if c in df.columns]
    df = df[keep_cols].copy()

    # 类型转换
    df["CHR"] = pd.to_numeric(df["CHR"], errors="coerce")
    df["BP"]  = pd.to_numeric(df["BP"], errors="coerce")
    df["P"]   = pd.to_numeric(df["P"], errors="coerce")
    if "MAF" in df.columns:
        df["MAF"] = pd.to_numeric(df["MAF"], errors="coerce")

    return df


def qc_one_trait(trait, fname):
    in_path = os.path.join(IN_DIR, fname)
    print(f"\n=== {trait}: {in_path} ===")

    df = load_and_standardize(in_path)
    n0 = len(df)
    print("原始 SNP 数:", n0)

    # 1) 去除缺失关键字段
    df = df.dropna(subset=["CHR","BP","P"])
    # 2) 只保留 1–20 号染色体
    df = df[df["CHR"].between(1, 20)]
    # 3) 去掉 P<=0 或 P>1 的
    df = df[(df["P"] > 0) & (df["P"] <= 1)]

    # 4) MAF 过滤（如果有这一列）
    if "MAF" in df.columns:
        before = len(df)
        df = df[(df["MAF"] >= MAF_MIN) & (df["MAF"] <= 0.5)]
        print(f"MAF 过滤后 SNP 数: {len(df)} (去掉 {before - len(df)})")
    else:
        print("注意：当前文件没有 MAF 列，暂时跳过 MAF 过滤。")

    # 5) 去除重复 SNP（按 chr:bp:A1:A2 组合去重）
    df["SNP_ID"] = df["SNP"].astype(str)
    key_cols = ["CHR","BP","A1","A2"]
    df["KEY"] = df[key_cols].astype(str).agg(":".join, axis=1)
    before = len(df)
    df = df.drop_duplicates(subset=["KEY"])
    print(f"按 chr:bp:A1:A2 去重后 SNP 数: {len(df)} (去掉 {before - len(df)})")

    # 输出干净版
    out_path = os.path.join(OUT_DIR, f"{trait}.qc.tsv")
    df.to_csv(out_path, sep="\t", index=False)
    print("QC 后文件已保存:", out_path)

    # 顺便输出一个简单的 lambda（用中位数近似）
    from numpy import median
    # χ^2 统计量
    chisq = (df["BETA"] / df["SE"])**2
    lam = chisq.median() / 0.456  # 0.456 是 1df χ^2 分布的中位数
    print(f"近似 lambda_GC: {lam:.3f}")

    return lam, len(df)


if __name__ == "__main__":
    summary = []
    for trait, fname in TRAITS.items():
        lam, n = qc_one_trait(trait, fname)
        summary.append((trait, n, lam))

    print("\n=== Wild GWAS QC 概览 ===")
    for trait, n, lam in summary:
        print(f"{trait}: {n} SNP, lambda≈{lam:.3f}")
