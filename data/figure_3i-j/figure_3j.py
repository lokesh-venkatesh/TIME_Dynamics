# plot_ifngamma.py

import pandas as pd
import matplotlib.pyplot as plt

# Load data
df = pd.read_csv("data/figure_3i-j/modified_run_oscillations.csv")

# Extract variables
S = df["S"]
IFNgamma = df["IFNgamma"]

# Create plot
plt.figure(figsize=(6,5))
plt.plot(IFNgamma, S, linewidth=2)

# Labels
plt.xlabel("IFN-γ")
plt.ylabel("S")

# Optional formatting similar to your figure
plt.ticklabel_format(style='sci', axis='x', scilimits=(6,6))

plt.grid(True)
plt.tight_layout()

# Show plot
plt.show()