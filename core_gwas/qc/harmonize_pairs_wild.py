# Generalized from soybean project-specific path layout.
#!/usr/bin/env python3
import os
import glob
import pandas as pd

WORKDIR = "data/input/workflow"
os.chdir(WORKDIR)

QC_DIR = "05_qc_wild"
CLUMP_DIR = "07_clumped_wild"
OUT_DIR = "08_harmonized_wild"
os.makedirs(OUT_DIR, exist_ok=True)

TRAITS = ["100SW", "Protein", "Oil", "log_ratio"]


def load_qc(trait):
    """读入 QC 后 summary."""
    path = os.path.join(QC_DIR, f"{trait}.qc.tsv")
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    df = pd.read_csv(path, sep="\t")
    needed = ["SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P"]
    missing = [c for c in needed if c not in df.columns]
    if missing:
        raise ValueError(f"{path} 缺少列: {missing}")
    return df


def find_clump_file(trait):
    """
    在 07_clumped_wild 里自动查找 {trait}.P*.clumped。
    如果有多个，优先选择 P 最小的那个（阈值最严格）。
    """
    pattern = os.path.join(CLUMP_DIR, f"{trait}.P*.clumped")
    files = glob.glob(pattern)
    if not files:
        print(f"  没找到 {pattern}，该 trait 暂时没有 clumped 结果。")
        return None

    # 尝试解析中间的 P 值，选最小的
    def parse_p(f):
        # 文件名形如 100SW.P0.0001.clumped 或 100SW.P1e-04.clumped
        base = os.path.basename(f)
        mid = base.split(".P", 1)[1].rsplit(".clumped", 1)[0]
        try:
            return float(mid)
        except Exception:
            # 解析失败就返回 1，当成很宽松
            return 1.0

    files_sorted = sorted(files, key=parse_p)
    chosen = files_sorted[0]
    print(f"  使用 clumped 文件: {os.path.basename(chosen)}")
    return chosen


def load_clumped_lead_snps(trait):
    """读入 clumped 文件中的 lead SNP 列表."""
    path = find_clump_file(trait)
    if path is None:
        return []

    df = pd.read_csv(path, delim_whitespace=True, comment="#")
    if "SNP" not in df.columns:
        raise ValueError(f"{path} 中没有 SNP 列")
    snps = df["SNP"].dropna().astype(str).unique().tolist()
    print(f"  {trait} clumped lead SNP 数: {len(snps)}")
    return snps


def is_palindromic(a1, a2):
    pair = (str(a1) + str(a2)).upper()
    return pair in ("AT", "TA", "CG", "GC")


def harmonize_pair(exposure, outcome):
    """
    对一对暴露/结局进行等位基因对齐，生成 MR 用 harmonized 表。
    输出：08_harmonized_wild/{exposure}_to_{outcome}.harmonized.tsv
    """
    print(f"\n=== {exposure} -> {outcome} ===")
    lead_snps = load_clumped_lead_snps(exposure)
    if len(lead_snps) == 0:
        print(f"  暴露 {exposure} 没有 clumped 工具 SNP，跳过。")
        return

    exp_df = load_qc(exposure)
    out_df = load_qc(outcome)

    # 1. 暴露端只保留 lead SNP
    exp_sub = exp_df[exp_df["SNP"].isin(lead_snps)].copy()
    print(f"  在 {exposure} GWAS 中找到的工具 SNP: {len(exp_sub)}")

    # 2. 与结局 GWAS 合并
    merged = exp_sub.merge(out_df, on="SNP", suffixes=("_exp", "_out"), how="inner")
    print(f"  在 {outcome} GWAS 中也存在的工具 SNP: {len(merged)}")
    if merged.empty:
        print("  暴露/结局没有共同 SNP，跳过。")
        return

    # 3. 丢掉 chr/pos 不一致的 SNP
    same_pos = (merged["CHR_exp"] == merged["CHR_out"]) & (merged["BP_exp"] == merged["BP_out"])
    if not same_pos.all():
        n_mis = (~same_pos).sum()
        print(f"  有 {n_mis} 个 SNP 在 chr/pos 上不一致，将丢弃。")
        merged = merged[same_pos].copy()

    # 4. 对齐等位基因（same / flip）
    a1e, a2e = merged["A1_exp"], merged["A2_exp"]
    a1o, a2o = merged["A1_out"], merged["A2_out"]

    same = (a1e == a1o) & (a2e == a2o)
    flip = (a1e == a2o) & (a2e == a1o)
    both = same | flip

    if not both.all():
        n_bad = (~both).sum()
        print(f"  有 {n_bad} 个 SNP A1/A2 无法简单对齐，丢弃。")
        merged = merged[both].copy()
        a1e, a2e = merged["A1_exp"], merged["A2_exp"]
        a1o, a2o = merged["A1_out"], merged["A2_out"]
        same = (a1e == a1o) & (a2e == a2o)
        flip = (a1e == a2o) & (a2e == a1o)

    beta_out = merged["BETA_out"].copy()
    se_out = merged["SE_out"].copy()
    beta_out[flip] = -beta_out[flip]

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

    # 5. 剔除 palindromic SNP（A/T, T/A, C/G, G/C）
    palin_mask = harm.apply(lambda r: is_palindromic(r["A1"], r["A2"]), axis=1)
    n_palin = palin_mask.sum()
    if n_palin > 0:
        print(f"  剔除 palindromic SNP 数: {n_palin}")
        harm = harm[~palin_mask].copy()

    n_final = len(harm)
    print(f"  最终可用工具 SNP 数: {n_final}")
    if n_final == 0:
        print("  无可用 SNP，未写出文件。")
        return

    out_path = os.path.join(OUT_DIR, f"{exposure}_to_{outcome}.harmonized.tsv")
    harm.to_csv(out_path, sep="\t", index=False)
    print("  已写出:", out_path)


if __name__ == "__main__":
    # 你关心的几对暴露/结局
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
