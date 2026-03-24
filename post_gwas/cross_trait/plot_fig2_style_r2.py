import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# ===============================
# Global style (Nature-like)
# ===============================
sns.set_style("white")
sns.set_context(
    "paper",
    rc={
        "axes.linewidth": 1.2,
        "xtick.major.width": 1.2,
        "ytick.major.width": 1.2
    }
)

# ===============================
# Load data
# ===============================
dat = pd.read_csv("MVMR_input_matrix.txt", sep="\t")

# ---------- helper ----------
def r2(x, y):
    return np.corrcoef(x, y)[0, 1] ** 2

def conditional_r2(y, x, z):
    # y ~ z, then correlate residuals with x
    coef = np.polyfit(z, y, 1)
    resid = y - (coef[0] * z + coef[1])
    return r2(resid, x)

# ===============================
# Fig.2b-like: size as outcome
# ===============================
df_size = pd.DataFrame({
    "Model": [
        "Oil",
        "Oil | Protein",
        "Protein",
        "Protein | Oil"
    ],
    "r2": [
        r2(dat["beta_size"], dat["beta_oil"]),
        conditional_r2(dat["beta_size"], dat["beta_oil"], dat["beta_protein"]),
        r2(dat["beta_size"], dat["beta_protein"]),
        conditional_r2(dat["beta_size"], dat["beta_protein"], dat["beta_oil"])
    ]
})

plt.figure(figsize=(5.5, 4))
ax = sns.stripplot(
    data=df_size,
    x="Model",
    y="r2",
    jitter=0.05,
    size=7,
    color="black"
)

ax.set_ylabel(r"$r^2$ (variance explained in size)")
ax.set_xlabel("")
ax.set_ylim(0, 1)
ax.set_title("Marginal and conditional variance explained in size")

plt.xticks(rotation=20)
sns.despine()
plt.tight_layout()
plt.savefig("Fig2b_style_size_clean.png", dpi=300)
plt.close()

# ===============================
# Fig.2c-like: oil / protein as outcome
# ===============================
df_trait = pd.DataFrame({
    "Outcome": (
        ["Oil"] * 4 +
        ["Protein"] * 4
    ),
    "Model": [
        "Size",
        "Protein",
        "Size | Protein",
        "Protein | Size",
        "Size",
        "Oil",
        "Size | Oil",
        "Oil | Size"
    ],
    "r2": [
        # oil
        r2(dat["beta_oil"], dat["beta_size"]),
        r2(dat["beta_oil"], dat["beta_protein"]),
        conditional_r2(dat["beta_oil"], dat["beta_size"], dat["beta_protein"]),
        conditional_r2(dat["beta_oil"], dat["beta_protein"], dat["beta_size"]),
        # protein
        r2(dat["beta_protein"], dat["beta_size"]),
        r2(dat["beta_protein"], dat["beta_oil"]),
        conditional_r2(dat["beta_protein"], dat["beta_size"], dat["beta_oil"]),
        conditional_r2(dat["beta_protein"], dat["beta_oil"], dat["beta_size"])
    ]
})

g = sns.catplot(
    data=df_trait,
    x="Model",
    y="r2",
    col="Outcome",
    kind="strip",
    jitter=0.05,
    color="black",
    height=4,
    aspect=1.1,
    sharey=True
)

g.set_axis_labels("", r"$r^2$ (variance explained)")
g.set(ylim=(0, 1))

for ax in g.axes.flat:
    for label in ax.get_xticklabels():
        label.set_rotation(25)

g.fig.suptitle(
    "Marginal and conditional variance explained in oil and protein",
    y=1.05
)

sns.despine()
plt.tight_layout()
plt.savefig("Fig2c_style_oil_protein_clean.png", dpi=300)
plt.close()

print("✓ Saved:")
print("  Fig2b_style_size_clean.png")
print("  Fig2c_style_oil_protein_clean.png")
