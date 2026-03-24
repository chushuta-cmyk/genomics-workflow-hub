# Generalized from soybean project-specific path layout.
#!/usr/bin/env python3
import os
import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

WORKDIR = "data/input/workflow"
os.chdir(WORKDIR)

HARM_DIR = "08_harmonized_wild"
RES_PATH = "11_mr_results_wild/mr_ivw_egger_wild_results.tsv"
PLOT_DIR = "12_mr_plots_wild"
os.makedirs(PLOT_DIR, exist_ok=True)

# 重点画图的组合
PLOT_PAIRS = [
    ("100SW", "Oil"),
    ("100SW", "Protein"),
    ("Oil", "Protein"),
    ("Protein", "Oil"),
]

# 读 MR 汇总结果
res = pd.read_csv(RES_PATH, sep="\t")

def get_mr_estimates(exposure, outcome):
    row = res[(res["exposure"] == exposure) & (res["outcome"] == outcome)]
    if row.empty:
        return None
    row = row.iloc[0]
    return {
        "beta_ivw": row["beta_ivw"],
        "se_ivw": row["se_ivw"],
        "beta_egger": row["beta_egger"],
        "se_egger": row["se_egger"],
        "alpha_egger": row["alpha_egger"],
        "se_alpha": row["se_alpha"],
        "n_snp": int(row["n_snp"]),
        "p_ivw": row["p_ivw"],
        "p_egger": row["p_egger"],
        "p_alpha": row["p_alpha"],
    }

def plot_scatter(exposure, outcome):
    fname = f"{exposure}_to_{outcome}.harmonized.tsv"
    path = os.path.join(HARM_DIR, fname)
    if not os.path.exists(path):
        print("找不到", path)
        return

    est = get_mr_estimates(exposure, outcome)
    if est is None:
        print("没有 MR 结果:", exposure, "->", outcome)
        return

    df = pd.read_csv(path, sep="\t")
    bx = df["BETA_exp"].values
    by = df["BETA_out"].values

    plt.figure(figsize=(5, 5))
    plt.scatter(bx, by, alpha=0.4, edgecolor="none")

    # x 轴范围
    x_min = bx.min()
    x_max = bx.max()
    x_pad = 0.05 * (x_max - x_min)
    xs = np.linspace(x_min - x_pad, x_max + x_pad, 100)

    # IVW 线（强制过原点）
    beta_ivw = est["beta_ivw"]
    ys_ivw = beta_ivw * xs
    plt.plot(xs, ys_ivw, color="red", label=f"IVW β={beta_ivw:.3f}")

    # Egger 线（带截距）
    alpha_eg = est["alpha_egger"]
    beta_eg = est["beta_egger"]
    ys_eg = alpha_eg + beta_eg * xs
    plt.plot(xs, ys_eg, color="blue", linestyle="--",
             label=f"Egger β={beta_eg:.3f}, α={alpha_eg:.3f}")

    plt.axhline(0, color="grey", linewidth=0.5)
    plt.axvline(0, color="grey", linewidth=0.5)

    plt.xlabel(f"{exposure} SNP effect (BETA_exp)")
    plt.ylabel(f"{outcome} SNP effect (BETA_out)")
    plt.title(f"{exposure} -> {outcome} (N SNP={est['n_snp']})")
    plt.legend()

    out_png = os.path.join(PLOT_DIR, f"{exposure}_to_{outcome}_scatter.png")
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()
    print("已保存图像:", out_png)


if __name__ == "__main__":
    print("=== 画 MR 散点图 + 回归线（wild 组） ===")
    for ex, out in PLOT_PAIRS:
        print(f"绘制 {ex} -> {out}")
        plot_scatter(ex, out)
