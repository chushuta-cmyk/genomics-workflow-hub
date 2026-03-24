# Generalized from soybean project-specific path layout.
import os
import glob
import pandas as pd



#输入目录和输出目录
WORKDIR = "data/input/workflow"
os.chdir(WORKDIR)

QC_DIR = "05_qc_wild"
CLUMP_DIR = "07_clumped_wild"
OUT_DIR = "08_harmonized_wild"
os.makedirs(OUT_DIR, exist_ok=True)

TRAITS = ["100SW", "Protein", "Oil", "log_ratio"]

#读 QC 后 GWAS 的函数
def load_qc(trait):
    path = os.path.join(QC_DIR, "{0}.qc.tsv".format(trait))
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    df = pd.read_csv(path, sep="\t")
    needed = ["SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P"]
    missing = [c for c in needed if c not in df.columns]
    if missing:
        raise ValueError("{0} 缺少列: {1}".format(path, missing))
    return df

#找到每个暴露的 clumped 文件并读出 lead SNP
def find_clump_file(trait):
    pattern = os.path.join(CLUMP_DIR, "{0}.P*.clumped".format(trait))
    files = glob.glob(pattern)
    if not files:
        print("  没找到 clumped 文件模式:", pattern)
        return None

    def parse_p(fname):
        base = os.path.basename(fname)
        mid = base.split(".P", 1)[1].rsplit(".clumped", 1)[0]
        try:
            return float(mid)
        except Exception:
            return 1.0

    files_sorted = sorted(files, key=parse_p)
    chosen = files_sorted[0]
    print("  使用 clumped 文件:", os.path.basename(chosen))
    return chosen

def load_lead_snps(trait):
    path = find_clump_file(trait)
    if path is None:
        return []

    df = pd.read_csv(path, delim_whitespace=True, comment="#")
    if "SNP" not in df.columns:
        raise ValueError("{0} 中没有 SNP 列".format(path))
    snps = df["SNP"].dropna().astype(str).unique().tolist()
    print("  {0}: clumped lead SNP 数 = {1}".format(trait, len(snps)))
    return snps

#判断 palindromic SNP 的小函数
def is_palindromic(a1, a2):
    pair = (str(a1) + str(a2)).upper()
    return pair in ("AT", "TA", "CG", "GC")


#核心函数：对齐一对暴露/结局
def harmonize_pair(exposure, outcome):
    print("\n=== {0} -> {1} ===".format(exposure, outcome))
    lead = load_lead_snps(exposure)
    if len(lead) == 0:
        print("  暴露 {0} 没有 clumped 工具 SNP，跳过。".format(exposure))
        return

    exp_df = load_qc(exposure)
    out_df = load_qc(outcome)

#用暴露工具在两个 GWAS 中取交集
    exp_sub = exp_df[exp_df["SNP"].isin(lead)].copy()
    print("  在 {0} GWAS 中找到的工具 SNP: {1}".format(exposure, len(exp_sub)))

    merged = exp_sub.merge(out_df, on="SNP", suffixes=("_exp", "_out"), how="inner")
    print("  在 {1} GWAS 中也存在的工具 SNP: {0}".format(len(merged), outcome))
    if merged.empty:
        print("  暴露/结局没有交集 SNP，跳过。")
        return

#确保坐标一致
    same_pos = (merged["CHR_exp"] == merged["CHR_out"]) & (merged["BP_exp"] == merged["BP_out"])
    if not same_pos.all():
        n_mis = (~same_pos).sum()
        print("  {0} 个 SNP 在 chr/pos 上不一致，将丢弃。".format(n_mis))
        merged = merged[same_pos].copy()

#等位基因对齐：same / flip
    a1e, a2e = merged["A1_exp"], merged["A2_exp"]
    a1o, a2o = merged["A1_out"], merged["A2_out"]

    same = (a1e == a1o) & (a2e == a2o)
    flip = (a1e == a2o) & (a2e == a1o)
    both = same | flip

    if not both.all():
        n_bad = (~both).sum()
        print("  {0} 个 SNP A1/A2 无法简单对齐，将丢弃。".format(n_bad))
        merged = merged[both].copy()
        a1e, a2e = merged["A1_exp"], merged["A2_exp"]
        a1o, a2o = merged["A1_out"], merged["A2_out"]
        same = (a1e == a1o) & (a2e == a2o)
        flip = (a1e == a2o) & (a2e == a1o)

    beta_out = merged["BETA_out"].copy()
    se_out = merged["SE_out"].copy()
    beta_out[flip] = -beta_out[flip]

#组装 harmonized 表
    harm = pd.DataFrame({
        "SNP": merged["SNP"],
        "CHR": merged["CHR_exp"],
        "BP": merged["BP_exp"],
        "A1": merged["A1_exp"],
        "A2": merged["A2_exp"],
        "BETA_exp": merged["BETA_exp"],
        "SE_exp": merged["SE_exp"],
        "BETA_out": beta_out,
        "SE_out": se_out,
    })

#去掉 palindromic SNP 和缺失
    palin_mask = harm.apply(lambda r: is_palindromic(r["A1"], r["A2"]), axis=1)
    n_palin = int(palin_mask.sum())
    if n_palin > 0:
        print("  剔除 palindromic SNP 数:", n_palin)
        harm = harm[~palin_mask].copy()

    before = len(harm)
    harm = harm.dropna(subset=["BETA_exp", "SE_exp", "BETA_out", "SE_out"])
    n_drop = before - len(harm)
    if n_drop > 0:
        print("  因缺失 beta/SE 剔除 SNP 数:", n_drop)


#写出结果
    n_final = len(harm)
    print("  最终可用工具 SNP 数:", n_final)
    if n_final == 0:
        print("  无可用 SNP，不写出文件。")
        return

    out_path = os.path.join(OUT_DIR, "{0}_to_{1}.harmonized.tsv".format(exposure, outcome))
    harm.to_csv(out_path, sep="\t", index=False)
    print("  已写出:", out_path)

#主程序
if __name__ == "__main__":
    pairs = [
        ("100SW", "Oil"),
        ("Oil", "100SW"),
        ("100SW", "Protein"),
        ("Protein", "100SW"),
        ("100SW", "log_ratio"),
        ("log_ratio", "100SW"),
        ("Oil", "Protein"),
        ("Protein", "Oil"),
    ]

    for ex, out in pairs:
        harmonize_pair(ex, out)
