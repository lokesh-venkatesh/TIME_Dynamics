import pandas as pd
import matplotlib.pyplot as plt

# ===============================
# 1. Import CSV file
# ===============================
file_path = "simulation_results.csv"  # Replace with your file name
df = pd.read_csv(file_path)

time = df["time"]

variables = [
    "S", "SR", "C", "CR",
    "M1", "M2",
    "TH1", "TH2",
    "TC", "Treg",
    "IL10", "IFNgamma", "IL2"
]

# ===============================
# 2. Create 13 separate plots
# ===============================
for var in variables:
    
    plt.figure(figsize=(8, 5))
    
    plt.plot(
        time,
        df[var],
        linewidth=2,
        color="tab:blue"
    )
    
    plt.xlabel("Time", fontsize=11)
    plt.ylabel(var, fontsize=11)
    plt.title(f"Time Series of {var}", fontsize=13)
    plt.grid(True, linestyle="--", alpha=0.6)
    
    # Optional: Use log scale for better visibility
    # Uncomment if needed
    # plt.yscale("log")
    
    plt.tight_layout()
    plt.savefig(f"initial_results/Variable {var} plotted versus time.png", dpi=300)
    plt.close()
