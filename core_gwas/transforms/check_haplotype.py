# Generalized from soybean project-specific path layout.
import pandas as pd
from pathlib import Path

f = Path("results/cultivated/fine_mapping/log_ratio_17_37608516/haplotype_groups.tsv")
df = pd.read_csv(f, sep="\t")

print(df.head())
print(df.columns.tolist())
print(df.dtypes)
print(df["group"].value_counts(dropna=False))

ycol = "log_Oil_Prot_ratio"
print("Non-null phenotype:", df[ycol].notna().sum())
print("Unique phenotype values:", df[ycol].nunique(dropna=True))
print(df.groupby("group")[ycol].describe())