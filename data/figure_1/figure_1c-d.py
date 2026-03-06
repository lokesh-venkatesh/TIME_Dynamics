import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

# Load data
df = pd.read_csv("./data/figure_1/default_run_results.csv")

plt.figure(figsize=(7,5))

"""
plt.plot(df["time"], df["S"], label="S", linewidth=2)
plt.plot(df["time"], df["SR"], label="SR", linewidth=2)
plt.plot(df["time"], df["C"], label="C", linewidth=2)
plt.plot(df["time"], df["CR"], label="CR", linewidth=2)
plt.plot(df["time"], df["M1"], label="M1", linewidth=2)
plt.plot(df["time"], df["M2"], label="M2", linewidth=2)
plt.plot(df["time"], df["TH1"], label="TH1", linewidth=2)
plt.plot(df["time"], df["TH2"], label="TH2", linewidth=2)
plt.plot(df["time"], df["TC"], label="TC", linewidth=2)
plt.plot(df["time"], df["Treg"], label="Treg", linewidth=2)
"""
plt.plot(df["time"], df["IL10"], color="red", linewidth=2, label="IL10")
plt.plot(df["time"], df["IFNgamma"], color="green", linewidth=2, label="IFN-γ")
plt.plot(df["time"], df["IL2"], color="black", linewidth=2, label="IL2")

plt.xlabel("Time (days)", fontsize=12)
plt.ylabel("Cell Density (cells/ml)", fontsize=12)
plt.legend(ncol=2, fontsize=8)
plt.xlim(left=0)
plt.ylim(bottom=0)

plt.tight_layout()
# plt.savefig("figure1/figure_1c.png", dpi=300)
# plt.close()
plt.show()