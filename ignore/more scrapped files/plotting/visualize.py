
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

def plot_dynamics(file_path, title, output_name):
    if not os.path.exists(file_path):
        print(f"File {file_path} not found. Run Julia scripts first.")
        return
        
    df = pd.read_csv(file_path)
    
    plt.figure(figsize=(12, 8))
    
    # Plotting Tumor Sub-populations
    plt.subplot(2, 1, 1)
    plt.plot(df['time'], df['S'], label='Stem (S)')
    plt.plot(df['time'], df['Sr'], label='Resistant Stem (Sr)')
    plt.plot(df['time'], df['C'], label='Cancer (C)')
    plt.plot(df['time'], df['Cr'], label='Resistant Cancer (Cr)')
    plt.yscale('log')
    plt.ylabel('Cell Count (Log Scale)')
    plt.title(f'{title} - Tumor Populations')
    plt.legend()
    plt.grid(True, alpha=0.3)

    # Plotting Cytokines
    plt.subplot(2, 1, 2)
    plt.plot(df['time'], df['IL10'], label='IL-10')
    plt.plot(df['time'], df['IFNg'], label='IFN-gamma')
    plt.plot(df['time'], df['IL2'], label='IL-2')
    plt.xlabel('Time (Days)')
    plt.ylabel('Concentration')
    plt.title(f'{title} - Cytokine Dynamics')
    plt.legend()
    plt.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_name)
    print(f"Figure saved as {output_name}")

if __name__ == "__main__":
    sns.set_theme(style="whitegrid")
    plot_dynamics("baseline_trajectory.csv", "Baseline Growth", "baseline_plot.png")
    plot_dynamics("treatment_protocol2_trajectory.csv", "Protocol 2 Treatment", "treatment_plot.png")
