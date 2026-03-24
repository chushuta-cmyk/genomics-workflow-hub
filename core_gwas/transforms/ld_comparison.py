# Generalized from soybean project-specific path layout.
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt

f = Path("results/comparison/ld_decay_comparison_summary.tsv")
df = pd.read_csv(f, sep="\t")

plt.figure(figsize=(7,5))
plt.plot(df["distance_bp"], df["wild_mean_r2"], marker="o", label="Wild")
plt.plot(df["distance_bp"], df["cultivated_mean_r2"], marker="s", label="Cultivated")
plt.xscale("log")
plt.xlabel("Distance (bp)")
plt.ylabel("Mean $r^2$")
plt.title("LD decay comparison (fixed distance points)")
plt.legend()
plt.tight_layout()
plt.savefig("results/comparison/ld_decay_comparison_fixed.png", dpi=300)
print(df)