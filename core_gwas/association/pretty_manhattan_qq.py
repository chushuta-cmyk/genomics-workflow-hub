#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# 兼容老版本 matplotlib
try:
    plt.style.use("seaborn")
except Exception:
    pass


def read_gemma_assoc(path):
    df = pd.read_csv(path, delim_whitespace=True, low_memory=False)
    df = df.rename(columns={
        "chr": "CHR",
        "rs": "SNP",
        "ps": "BP",
        "p_wald": "P"
    })
    df = df[df["P"] > 0].copy()
    return df


def manhattan_qq(
    assoc_file,
    title="Trait",
    suggestiveline=5e-5,
    signifline=5e-8,
    save_prefix="plots/out"
):
    df = read_gemma_assoc(assoc_file)

    # 1. 只保留 1–20 号染色体
    df["CHR"] = pd.to_numeric(df["CHR"], errors="coerce")
    df = df[df["CHR"].isin(range(1, 21))].copy()

    # 2. 按 1–20 排序
    df = df.sort_values(["CHR", "BP"])
    chr_groups = df.groupby("CHR")

    # 计算累计坐标
    chr_offsets = {}
    current_offset = 0
    ticks = []
    labels = []

    for chrom in range(1, 21):
        if chrom not in chr_groups.groups:
            continue
        group = chr_groups.get_group(chrom)
        chr_offsets[chrom] = current_offset
        max_bp = group["BP"].max()
        ticks.append(current_offset + max_bp / 2)
        labels.append(str(chrom))
        current_offset += max_bp

    df["BP_cum"] = df.apply(lambda r: r["BP"] + chr_offsets[r["CHR"]], axis=1)
    df["minus_log10_p"] = -np.log10(df["P"])

    # 画布
    fig = plt.figure(figsize=(10, 8))
    gs = GridSpec(2, 1, height_ratios=[2.5, 1.5])

    colors = ["#4C72B0", "#55A868"]

    # 曼哈顿图
    ax1 = fig.add_subplot(gs[0])
    for i, chrom in enumerate(range(1, 21)):
        if chrom not in chr_groups.groups:
            continue
        group = df[df["CHR"] == chrom]
        ax1.scatter(
            group["BP_cum"],
            group["minus_log10_p"],
            c=colors[i % 2],
            s=6,
            alpha=0.7,
            linewidths=0
        )

    if suggestiveline:
        ax1.axhline(-np.log10(suggestiveline), color="orange", linestyle="--", linewidth=1)
    if signifline:
        ax1.axhline(-np.log10(signifline), color="red", linestyle="-.", linewidth=1)

    ax1.set_xticks(ticks)
    ax1.set_xticklabels(labels, fontsize=8)
    ax1.set_ylabel(r"$-log_{10}(P)$")
    ax1.set_title(f"{title} (Manhattan)", fontsize=12)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    # QQ 图
    ax2 = fig.add_subplot(gs[1])
    p_values = np.sort(df["P"].values)
    n = len(p_values)

    expected = -np.log10((np.arange(1, n + 1) - 0.5) / n)
    observed = -np.log10(p_values)

    ax2.scatter(expected, observed, c="#4C72B0", s=8, alpha=0.6, edgecolor="none")
    max_val = max(expected.max(), observed.max())
    ax2.plot([0, max_val], [0, max_val], color="red", linestyle="--", linewidth=1)

    ax2.set_xlabel("Expected -log10(P)")
    ax2.set_ylabel("Observed -log10(P)")
    ax2.set_title(f"{title} (QQ)", fontsize=12)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    fig.subplots_adjust(hspace=0.3)

    import os
    out_dir = os.path.dirname(save_prefix)
    if out_dir and not os.path.isdir(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    plt.savefig(f"{save_prefix}_manhattan_qq.png", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    traits = [
        {
            "file": "04_gemma_wild/wild_trait_1.assoc.txt",
            "title": "Wild 100SW",
            "prefix": "plots/wild_100SW"
        },
        {
            "file": "04_gemma_wild/wild_trait_2.assoc.txt",
            "title": "Wild Protein",
            "prefix": "plots/wild_protein"
        },
        {
            "file": "04_gemma_wild/wild_trait_3.assoc.txt",
            "title": "Wild Oil",
            "prefix": "plots/wild_oil"
        },
        {
            "file": "04_gemma_wild/wild_trait_5.assoc.txt",
            "title": "Wild log(Oil/Protein)",
            "prefix": "plots/wild_log_ratio"
        },
    ]

    for t in traits:
        print("Plotting:", t["title"])
        manhattan_qq(
            assoc_file=t["file"],
            title=t["title"],
            save_prefix=t["prefix"]
        )
