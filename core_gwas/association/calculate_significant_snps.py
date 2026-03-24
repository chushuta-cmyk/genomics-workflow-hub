# Generalized from soybean project-specific path layout.
import pandas as pd
from pathlib import Path

base = Path("results/wild/post_gwas/significant_snps")

traits = ["100SW","Protein","Oil","log_ratio"]

rows=[]

for t in traits:

    f = base / f"sig_{t}.tsv"

    df = pd.read_csv(f, sep="\t")

    rows.append({
        "trait":t,
        "num_sig_snps":len(df)
    })

res = pd.DataFrame(rows)

print(res)

res.to_csv("significant_snp_counts_fixed.tsv", sep="\t", index=False)