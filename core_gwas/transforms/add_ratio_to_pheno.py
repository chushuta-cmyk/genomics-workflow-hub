# add_ratio_to_pheno.py
import pandas as pd
import numpy as np

pheno = pd.read_csv("../data/phenotype_original.tsv", sep="\t")

# 假设原来的列名为 'Oil' 和 'Protein'；如果不一样改这里
pheno["Oil_Prot_ratio"] = pheno["Oil"] / pheno["Protein"]
pheno["log_Oil_Prot_ratio"] = (pheno["Oil"] / pheno["Protein"]).apply(lambda x: np.log(x))

pheno.to_csv("../data/phenotype_with_ratio.tsv", sep="\t", index=False)
print("✓ phenotype_with_ratio.tsv 已生成")
