# Generalized from soybean project-specific path layout.
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

base = Path("data/projects/soybean_analysis")
wild_dir = base / "wild/post_gwas/significant_snps"
cult_dir = base / "cultivated/post_gwas"

trait_map = {
    "100SW": ("sig_100SW.tsv", "sig_100SW.tsv"),
    "Protein": ("sig_Protein.tsv", "sig_Protein.tsv"),
    "Oil": ("sig_Oil.tsv", "sig_Oil.tsv"),
    "log_ratio": ("sig_log_ratio.tsv", "sig_log_ratio.tsv"),
}

rows = []
for trait, (wfile, cfile) in trait_map.items():
    w = pd.read_csv(wild_dir / wfile, sep="\t")
    c = pd.read_csv(cult_dir / cfile, sep="\t")
    rows.append({
        "trait": trait,
        "wild_sig_snps": len(w),
        "cultivated_sig_snps": len(c),
        "total_unique": len(set(w["rs"]).union(set(c["rs"])))
    })

df = pd.DataFrame(rows)
out_tsv = base / "comparison/gwas_comparison_summary_fixed.tsv"
df.to_csv(out_tsv, sep="\t", index=False)

x = np.arange(len(df))
width = 0.35

plt.figure(figsize=(8,5))
plt.bar(x - width/2, df["wild_sig_snps"], width, label="Wild")
plt.bar(x + width/2, df["cultivated_sig_snps"], width, label="Cultivated")
plt.xticks(x, df["trait"])
plt.ylabel("Number of significant SNPs")
plt.title("Significant SNP count comparison (fixed)")
plt.legend()
plt.tight_layout()
plt.savefig(base / "comparison/significant_snp_count_comparison_fixed.png", dpi=300)
print(df)