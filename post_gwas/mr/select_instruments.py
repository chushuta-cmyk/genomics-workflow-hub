# Generalized from soybean project-specific path layout.
#!/usr/bin/env python3
import os
import subprocess
import pandas as pd

# 工作目录
WORKDIR = "data/input/workflow"
os.chdir(WORKDIR)

QC_DIR = "05_qc_wild"
CLUMP_DIR = "07_clumped_wild"
os.makedirs(CLUMP_DIR, exist_ok=True)

# wild 103 个样本，用 P<1e-4 作为主分析阈值
P_THRESH = 1e-4

TRAITS = ["100SW", "Protein", "Oil", "log_ratio"]

# !!! 使用刚刚生成的、只含 chr1-20 的 LD 参考面板
BFILE = "01_plink/05_final/soybean_wild_chr1_20"

print("=== Wild 工具 SNP 选择与 clumping (P < {0:g}) ===".format(P_THRESH))

for trait in TRAITS:
    print("\n=== 处理性状: {0} ===".format(trait))

    # 1. 统计 QC 文件中 P<阈值 的 SNP 数（仅记录）
    qc_path = os.path.join(QC_DIR, "{0}.qc.tsv".format(trait))
    if not os.path.exists(qc_path):
        print("  找不到 QC 文件:", qc_path)
        continue

    df = pd.read_csv(qc_path, sep="\t")
    n_sig = int((df["P"] < P_THRESH).sum())
    print("  P<{0:g} 候选 SNP (未 clump 前) = {1}".format(P_THRESH, n_sig))

    # 2. 在整个 QC 文件上 clump，由 --clump-p1 控制 index SNP 阈值
    clump_out = os.path.join(CLUMP_DIR, "{0}.P0.0001".format(trait))
    cmd = [
        "plink",
        "--bfile", BFILE,
        "--allow-extra-chr",          # 允许非常规染色体名，避免 scaffold 报错
        "--clump", qc_path,
        "--clump-field", "P",
        "--clump-p1", str(P_THRESH),
        "--clump-p2", "1.0",
        "--clump-r2", "0.1",
        "--clump-kb", "250",
        "--out", clump_out
    ]
    print("  运行 PLINK clumping ...")
    try:
        res = subprocess.run(cmd, check=True,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        print("  clumping 完成:", clump_out + ".clumped")
    except subprocess.CalledProcessError as e:
        print("  clumping 失败，PLINK stderr 如下：")
        print(e.stderr.decode("utf-8"))
        continue

    # 3. 统计独立工具数量
    clump_path = clump_out + ".clumped"
    if os.path.exists(clump_path):
        df_clump = pd.read_csv(clump_path, delim_whitespace=True, comment="#")
        n_inst = len(df_clump)  # 一般一行一个 lead SNP
        print("  {0}: {1} 个独立工具 SNP ✓".format(trait, n_inst))
    else:
        print("  {0}: clumped 文件未生成 ✗".format(trait))

print("\n=== 汇总 ===")
for trait in TRAITS:
    clump_path = os.path.join(CLUMP_DIR, "{0}.P0.0001.clumped".format(trait))
    if os.path.exists(clump_path):
        df_clump = pd.read_csv(clump_path, delim_whitespace=True, comment="#")
        print("{0}: {1} 个独立工具 SNP".format(trait, len(df_clump)))
    else:
        print("{0}: 无 clumped 文件".format(trait))
