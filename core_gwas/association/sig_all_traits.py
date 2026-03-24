# Generalized from soybean project-specific path layout.
import pandas as pd

df = pd.read_csv("results/wild/post_gwas/significant_snps/sig_all_traits.tsv", sep="\t")

df = df.sort_values(["chr","ps"])

window = 250000

loci=[]

for trait in df.trait.unique():

    sub=df[df.trait==trait]

    sub=sub.sort_values(["chr","ps"])

    current_chr=None
    current_end=0

    for _,r in sub.iterrows():

        if r.chr!=current_chr or r.ps>current_end:

            loci.append({
                "trait":trait,
                "chr":r.chr,
                "start":r.ps,
                "end":r.ps
            })

            current_chr=r.chr
            current_end=r.ps+window

        else:

            loci[-1]["end"]=r.ps

loci_df=pd.DataFrame(loci)

print(loci_df.groupby("trait").size())

loci_df.to_csv("trait_specific_loci.tsv",sep="\t",index=False)