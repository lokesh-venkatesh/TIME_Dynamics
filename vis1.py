import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# ===============================
# 1. Import CSV file
# ===============================
# Replace with your actual file name
file_path = "simulation_results.csv"

df = pd.read_csv(file_path)

# ===============================
# 2. Select time and 13 time series
# ===============================
time = df["time"]

variables = [
    "S", "SR", "C", "CR",
    "M1", "M2",
    "TH1", "TH2",
    "TC", "Treg",
    "IL10", "IFNgamma", "IL2"
]

# ===============================
# 3. Plot configuration
# ===============================
plt.figure(figsize=(14, 8))

# Generate distinct colors
colors = plt.cm.tab20(np.linspace(0, 1, len(variables)))

# Different line styles and markers
linestyles = ['-', '--', '-.', ':']
markers = ['o', 's', 'D', '^', 'v', '<', '>', 'p', '*', 'x', '+', 'h', '.']

# ===============================
# 4. Plot each time series
# ===============================
for i, var in enumerate(variables):
    plt.plot(
        time,
        df[var],
        label=var,
        color=colors[i],
        linestyle=linestyles[i % len(linestyles)],
        marker=markers[i % len(markers)],
        markevery=int(len(time)/20),  # reduce marker density
        linewidth=1,
        markersize=6
    )

# ===============================
# 5. Formatting
# ===============================
plt.xlabel("Time", fontsize=12)
plt.ylabel("Population / Concentration", fontsize=12)

# TO ENABLE OR NOT TO ENABLE
# plt.yscale("log")

plt.title("Simulation Results: 13 Coupled Time Series", fontsize=14)
plt.legend(loc="upper left", bbox_to_anchor=(1, 1))
plt.grid(True, linestyle="--", alpha=0.4)
plt.tight_layout()

# Optional: Use logarithmic y-scale if values span many orders
# plt.yscale("log")

plt.savefig("initial_results/Variables plotted together against time.png", dpi=300)
plt.close()