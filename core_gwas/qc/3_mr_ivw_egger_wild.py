# Generalized from soybean project-specific path layout.
#!/usr/bin/env python3
import os
import math
import numpy as np
import pandas as pd

WORKDIR = "data/input/workflow"
os.chdir(WORKDIR)

HARM_DIR = "08_harmonized_wild"
OUT_DIR = "11_mr_results_wild"
os.makedirs(OUT_DIR, exist_ok=True)

# 要分析的暴露->结局组合（和你 harmonize 一致）
PAIRS = [
    ("100SW", "Oil"),
    ("Oil", "100SW"),
    ("100SW", "Protein"),
    ("Protein", "100SW"),
    ("100SW", "log_ratio"),
    ("log_ratio", "100SW"),
    ("Oil", "Protein"),
    ("Protein", "Oil"),
]


def normal_pvalue(z):
    """双侧正态检验 p 值（不依赖 scipy）."""
    # 标准正态 CDF 近似: 0.5 * (1 + erf(z / sqrt(2)))
    p_one_side = 0.5 * (1.0 + math.erf(-abs(z) / math.sqrt(2.0)))
    return 2.0 * p_one_side


def ivw_mr(bx, by, se_by):
    """
    固定效应 IVW，无截距：
    回归 by ~ bx，权重 w = 1/se_by^2
    """
    w = 1.0 / (se_by ** 2)
    num = np.sum(w * bx * by)
    den = np.sum(w * bx * bx)
    if den <= 0:
        return np.nan, np.nan, np.nan

    beta_ivw = num / den
    se_ivw = math.sqrt(1.0 / den)
    z = beta_ivw / se_ivw
    p = normal_pvalue(z)
    return beta_ivw, se_ivw, p


def egger_mr(bx, by, se_by):
    """
    MR-Egger：加截距的加权线性回归 by = alpha + beta * bx + e
    权重 w = 1/se_by^2
    返回: beta, se_beta, p_beta, alpha, se_alpha, p_alpha
    """
    w = 1.0 / (se_by ** 2)
    X = np.vstack([np.ones_like(bx), bx]).T  # 列: [截距, bx]
    W = np.diag(w)
    XtW = X.T @ W
    XtWX = XtW @ X
    try:
        XtWX_inv = np.linalg.inv(XtWX)
    except np.linalg.LinAlgError:
        return [np.nan] * 6

    beta_hat = XtWX_inv @ (XtW @ by)
    alpha, beta = beta_hat[0], beta_hat[1]

    # 残差 & 残差方差
    y_hat = X @ beta_hat
    resid = by - y_hat
    # 自由度 n-2
    dof = max(len(by) - 2, 1)
    sigma2 = np.sum(w * resid * resid) / dof

    cov_beta = XtWX_inv * sigma2
    se_alpha = math.sqrt(cov_beta[0, 0])
    se_beta = math.sqrt(cov_beta[1, 1])

    z_alpha = alpha / se_alpha if se_alpha > 0 else np.nan
    z_beta = beta / se_beta if se_beta > 0 else np.nan
    p_alpha = normal_pvalue(z_alpha) if not np.isnan(z_alpha) else np.nan
    p_beta = normal_pvalue(z_beta) if not np.isnan(z_beta) else np.nan

    return beta, se_beta, p_beta, alpha, se_alpha, p_alpha


def run_mr_for_pair(exposure, outcome):
    fname = f"{exposure}_to_{outcome}.harmonized.tsv"
    path = os.path.join(HARM_DIR, fname)
    if not os.path.exists(path):
        print(f"  找不到 {path}，跳过。")
        return None

    df = pd.read_csv(path, sep="\t")
    n_snp = len(df)
    if n_snp < 3:
        print(f"  {fname}: SNP 数 ({n_snp}) 过少，跳过 MR。")
        return None

    bx = df["BETA_exp"].values.astype(float)
    by = df["BETA_out"].values.astype(float)
    se_by = df["SE_out"].values.astype(float)

    # 去掉 SE 为 0 或缺失的点
    mask = (~np.isnan(bx)) & (~np.isnan(by)) & (~np.isnan(se_by)) & (se_by > 0)
    bx = bx[mask]
    by = by[mask]
    se_by = se_by[mask]
    n_use = len(bx)
    if n_use < 3:
        print(f"  {fname}: 有效 SNP 数 ({n_use}) 过少，跳过 MR。")
        return None

    print(f"  {fname}: 用 {n_use} 个 SNP 做 MR 分析。")

    # IVW MR
    beta_ivw, se_ivw, p_ivw = ivw_mr(bx, by, se_by)

    # MR-Egger
    (beta_egger, se_egger, p_egger,
     alpha_egger, se_alpha, p_alpha) = egger_mr(bx, by, se_by)

    res = {
        "exposure": exposure,
        "outcome": outcome,
        "n_snp": n_use,
        "beta_ivw": beta_ivw,
        "se_ivw": se_ivw,
        "p_ivw": p_ivw,
        "beta_egger": beta_egger,
        "se_egger": se_egger,
        "p_egger": p_egger,
        "alpha_egger": alpha_egger,
        "se_alpha": se_alpha,
        "p_alpha": p_alpha,
    }
    return res


if __name__ == "__main__":
    results = []

    print("=== Wild 组 MR (IVW + MR-Egger) 分析 ===")
    for ex, out in PAIRS:
        print(f"\n>>> {ex} -> {out}")
        res = run_mr_for_pair(ex, out)
        if res is not None:
            results.append(res)

    if results:
        df_res = pd.DataFrame(results)
        out_csv = os.path.join(OUT_DIR, "mr_ivw_egger_wild_results.tsv")
        df_res.to_csv(out_csv, sep="\t", index=False)
        print("\n=== 结果汇总已写出 ===")
        print("  ", out_csv)
        print(df_res)
    else:
        print("\n没有可用结果。")
