import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

# Load data
df = pd.read_csv("./data/figure_1/default_run_results.csv")
last = df.iloc[-1]

CD4 = last["TH1"] + last["TH2"]
CD8 = last["TC"]
Treg = last["Treg"]

ratios = [
    CD4 / CD8,
    CD4 / Treg,
    CD8 / Treg
]

labels = ["CD4/CD8", "CD4/Treg", "CD8/Treg"]

plt.figure(figsize=(6,5))
plt.bar(labels, ratios, color="maroon")

plt.ylabel("Cell Ratio", fontsize=12)
plt.ylim(bottom=0)

plt.tight_layout()
# plt.savefig("figure1/figure_1e.png", dpi=300)
# plt.close()
plt.show()