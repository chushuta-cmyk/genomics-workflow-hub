import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

base = Path("/data03/karama/projects/soybean_analysis/comparison")
df = pd.read_csv(base / "shared_loci_summary_fixed.tsv", sep="\t")

x = np.arange(len(df))
width = 0.25

plt.figure(figsize=(10,6))
plt.bar(x - width, df["wild_total"], width, label="Wild total")
plt.bar(x, df["cultivated_total"], width, label="Cultivated total")
plt.bar(x + width, df["shared"], width, label="Shared")
plt.xticks(x, df["trait"])
plt.ylabel("Number of loci")
plt.title("Trait-specific locus comparison (fixed)")
plt.legend()
plt.tight_layout()
plt.savefig(base / "locus_count_comparison_fixed.png", dpi=300)