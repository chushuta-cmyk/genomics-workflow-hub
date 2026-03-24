import pandas as pd
from pathlib import Path

base = Path("/data03/karama/projects/soybean_analysis")

wild_dir = base / "wild/post_gwas"
cult_dir = base / "cultivated/post_gwas"

traits = ["100SW", "Protein", "Oil", "log_ratio"]

rows = []

WINDOW = 250000

for trait in traits:
    wild = pd.read_csv(wild_dir / f"locus_summary_{trait}.tsv", sep="\t")
    cult = pd.read_csv(cult_dir / f"locus_summary_{trait}.tsv", sep="\t")

    shared = 0
    matched = set()

    for i, w in wild.iterrows():
        for j, c in cult.iterrows():

            if j in matched:
                continue

            if str(w["chr"]) != str(c["chr"]):
                continue

            if abs(int(w["lead_snp_pos"]) - int(c["lead_snp_pos"])) <= WINDOW:
                shared += 1
                matched.add(j)
                break

    rows.append({
        "trait": trait,
        "wild_total": len(wild),
        "cultivated_total": len(cult),
        "shared": shared,
        "wild_unique": len(wild) - shared,
        "cultivated_unique": len(cult) - shared
    })

df = pd.DataFrame(rows)

print(df)

df.to_csv(base / "comparison/shared_loci_summary_fixed.tsv",
          sep="\t", index=False)