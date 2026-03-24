import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

base = Path("/data03/karama/projects/soybean_analysis")

wild_dir = base / "wild/post_gwas"
cult_dir = base / "cultivated/post_gwas"
shared_file = base / "comparison/shared_loci_summary.tsv"

traits = ["100SW", "Protein", "Oil", "log_ratio"]

rows = []
shared_df = pd.read_csv(shared_file, sep="\t")

for trait in traits:
    w = pd.read_csv(wild_dir / f"locus_summary_{trait}.tsv", sep="\t")
    c = pd.read_csv(cult_dir / f"locus_summary_{trait}.tsv", sep="\t")
    sub = shared_df[(shared_df["trait"] == trait) & (shared_df["type"] == "shared_pair")]
    rows.append({
        "trait": trait,
        "wild_total": len(w),
        "cultivated_total": len(c),
        "shared": len(sub),
        "wild_unique": len(w) - len(sub),
        "cultivated_unique": len(c) - len(sub),
    })

df = pd.DataFrame(rows)
df.to_csv(base / "comparison/locus_count_comparison_fixed.tsv", sep="\t", index=False)

x = np.arange(len(df))
width = 0.25

plt.figure(figsize=(10,6))
plt.bar(x - width, df["wild_total"], width, label="Wild total")
plt.bar(x, df["cultivated_total"], width, label="Cultivated total")
plt.bar(x + width, df["shared"], width, label="Shared")
plt.xticks(x, df["trait"])
plt.ylabel("Number of loci")
plt.title("Trait-specific locus count comparison (fixed)")
plt.legend()
plt.tight_layout()
plt.savefig(base / "comparison/locus_count_comparison_fixed.png", dpi=300)
print(df)