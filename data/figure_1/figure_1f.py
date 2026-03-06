import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

# Load data
df = pd.read_csv("./data/figure_1/default_run_results.csv")

initial = df.iloc[0]
final = df.iloc[-1]

def compute_ratios(row):
    CD4 = row["TH1"] + row["TH2"]
    CD8 = row["TC"]
    Treg = row["Treg"]
    return [
        CD4 / CD8,
        CD4 / Treg,
        CD8 / Treg,
        row["TH1"] / row["TH2"],
        row["M1"] / row["M2"]
    ]

ratios_initial = compute_ratios(initial)
ratios_final = compute_ratios(final)

labels = ["CD4/CD8", "CD4/Treg", "CD8/Treg", "TH1/TH2", "M1/M2"]

x = np.arange(len(labels))
width = 0.35

fig, ax = plt.subplots(figsize=(8,5))

# Better color palette
color_initial = "#4C72B0"  # blue
color_final = "#DD8452"    # orange

bars1 = ax.bar(x - width/2, ratios_initial, width,
               label="Disease free", color=color_initial)

bars2 = ax.bar(x + width/2, ratios_final, width,
               label="Cancer", color=color_final)

# Add value labels on top of bars
def add_labels(bars):
    for bar in bars:
        height = bar.get_height()
        ax.text(
            bar.get_x() + bar.get_width()/2,
            height,
            f"{height:.2f}",
            ha='center',
            va='bottom',
            fontsize=9
        )

add_labels(bars1)
add_labels(bars2)

ax.set_xticks(x)
ax.set_xticklabels(labels)

ax.set_ylabel("Cell Ratio", fontsize=12)
ax.set_title("Immune Cell Population Ratios: Disease-Free vs Cancer", fontsize=14)

ax.legend()
ax.set_ylim(bottom=0)

# Explicit biological definition note
fig.text(
    0.5,
    -0.02,
    "CD4 = TH1 + TH2; CD8 = Tc",
    ha="center",
    fontsize=10
)

plt.tight_layout()
# plt.savefig("figure1/figure_1f.png", dpi=300, bbox_inches="tight")
# plt.close()
plt.show()