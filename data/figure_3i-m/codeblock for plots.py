os.makedirs('data', exist_ok=True)
os.makedirs("data/figure_3i-m", exist_ok=True)
default_simulation_results = run_sweep({'muTh1Ck1': np.linspace(0, 100, 50)}, output_prefix='sweep_muTh1Ck1')

# Extract metrics from sweep results - iterate through dictionary items directly
param_values = []
th1_th2_ratios = []
fold_changes = {'S': [], 'SR': [], 'C': [], 'CR': []}

for param_val, df in sorted(default_simulation_results.items()): # Iterate through (key, value) pairs directly
    param_values.append(param_val)
    
    # Get final steady-state values
    final_th1 = df['TH1'].iloc[-1]
    final_th2 = df['TH2'].iloc[-1]
    
    # Calculate TH1/TH2 ratio
    th1_th2_ratio = final_th1 / final_th2 if final_th2 > 0 else np.nan
    th1_th2_ratios.append(th1_th2_ratio)
    
    # Calculate fold changes (final / initial) - handle zero initial conditions
    # print(df['S'].iloc[-1], df['SR'].iloc[-1], df['C'].iloc[-1], df['CR'].iloc[-1]) # Debug print for initial values
    # initial values = 23618154489.75116 21250.82882227321 750051321.053279 153613636.80424753
    initial_s, initial_sr, initial_c, initial_cr = 23618154489.75116, 21250.82882227321, 750051321.053279, 153613636.80424753

    # initial_s = df['S'].iloc[0]
    # initial_sr = df['SR'].iloc[0] if df['SR'].iloc[0] > 0 else 0 # np.nan # 1e-10
    # initial_c = df['C'].iloc[0] if df['C'].iloc[0] > 0 else 0 # np.nan # 1e-10
    # initial_cr = df['CR'].iloc[0] if df['CR'].iloc[0] > 0 else 0 # np.nan # 1e-10
    
    fold_changes['S'].append(df['S'].iloc[-1] / initial_s if initial_s > 0 else 1)
    fold_changes['SR'].append(df['SR'].iloc[-1] / initial_sr if initial_sr > 0 else 1)
    fold_changes['C'].append(df['C'].iloc[-1] / initial_c if initial_c > 0 else 1)
    fold_changes['CR'].append(df['CR'].iloc[-1] / initial_cr if initial_cr > 0 else 1) # else np.nan

print(f"Extracted {len(param_values)} parameter values")

# Plot 1: TH1/TH2 ratio
fig1, ax1 = plt.subplots(figsize=(8, 6))
ax1.plot(param_values, th1_th2_ratios, 'b-o', linewidth=2, markersize=6)
ax1.set_xlabel('μ$_{Th1Ck1}$', fontsize=12, fontweight='bold')
ax1.set_ylabel('TH1/TH2', fontsize=12, fontweight='bold')
# ax1.text(0.05, 0.95, 'l', transform=ax1.transAxes, fontsize=14, fontweight='bold', va='top')
ax1.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('data/figure_3l.png', dpi=300, bbox_inches='tight')
plt.show()
plt.close()

# Plot 2: Fold changes
fig2, ax2 = plt.subplots(figsize=(8, 6))
colors = {'S': 'blue', 'SR': 'green', 'C': 'red', 'CR': 'cyan'}
for species, color in colors.items():
    ax2.plot(param_values, fold_changes[species], color=color, linewidth=2, label=species, marker='o', markersize=4)

ax2.set_xlabel('μ$_{Th1Ck1}$', fontsize=12, fontweight='bold')
ax2.set_ylabel('Fold Change', fontsize=12, fontweight='bold')
ax2.text(0.05, 0.95, 'm', transform=ax2.transAxes, fontsize=14, fontweight='bold', va='top')
ax2.legend(loc='best', fontsize=10)
ax2.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('data/figure_3m.png', dpi=300, bbox_inches='tight')
plt.show()