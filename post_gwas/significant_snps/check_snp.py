import pandas as pd

df = pd.read_csv("sig_all_traits.tsv", sep="\t")

df = df.sort_values(["chr","ps"])

df = df.rename(columns={
    "ps":"pos",
    "p_wald":"p"
})

df.to_csv("all_significant_snps.tsv", sep="\t", index=False)

print("Total SNP:", len(df))
print(df.groupby("trait").size())