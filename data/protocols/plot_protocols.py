import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Load data
df_p1_no_resist = pd.read_csv('data/protocols/protocol_1_without_resistant.csv')
df_p1_resist = pd.read_csv('data/protocols/protocol_1_with_resistant.csv')
df_p2 = pd.read_csv('data/protocols/protocol_2_combined.csv')

# Define treatment phases for coloring
# Phase 1: 0-200 (Growth, black)
# Phase 2: 200-300 (Chemotherapy, green)
# Phase 3: 300-400 (Radiotherapy, red)
# Phase 4: 400+ (Treatment-free or maintenance, varies)

def get_phase_colors(time_array, protocol='p1'):
    """Determine color for each time point based on treatment phase"""
    colors = []
    for t in time_array:
        if t < 200:
            colors.append('black')
        elif t < 300:
            colors.append('green')
        elif t < 400:
            colors.append('red')
        elif protocol == 'p2' and t < 500:
            colors.append('magenta')
        else:
            colors.append('black')
    return colors

# Create figure with 4 rows x 3 columns
fig, axes = plt.subplots(4, 3, figsize=(16, 14))

# Define which column corresponds to which protocol and data
protocols_data = [
    (df_p1_no_resist, 'Protocol 1\nwithout resistant cells'),
    (df_p1_resist, 'Protocol 1\nwith resistant cells'),
    (df_p2, 'Protocol 2\nwith resistant cells')
]

# Define rows (species to plot)
species_info = [
    ('S', 'Stem Cells'),
    ('SR', 'Stem Resistant Cells'),
    ('C', 'Cancer Cells'),
    ('CR', 'Cancer Resistant Cells')
]

# Plot each subplot
for col_idx, (df, title) in enumerate(protocols_data):
    # Determine protocol type for coloring
    protocol_type = 'p2' if col_idx == 2 else 'p1'
    
    for row_idx, (species, species_title) in enumerate(species_info):
        ax = axes[row_idx, col_idx]
        
        # Get data
        time = df['time'].values
        y_data = df[species].values
        
        # Plot with color segments based on treatment phase
        phase_colors = get_phase_colors(time, protocol_type)
        
        # Plot as continuous line with varying colors
        for i in range(len(time) - 1):
            ax.plot(time[i:i+2], y_data[i:i+2], color=phase_colors[i], linewidth=2)
        
        # Set labels and title
        if row_idx == 0:
            ax.set_title(title, fontsize=12, fontweight='bold')
        if col_idx == 0:
            ax.set_ylabel(species_title, fontsize=11, fontweight='bold')
        
        ax.set_xlabel('Time (days)', fontsize=10)
        ax.grid(True, alpha=0.3)
        
        # Set x-axis limits based on protocol
        if protocol_type == 'p2':
            ax.set_xlim(0, 800)
        else:
            ax.set_xlim(0, 500)
        
        # Format y-axis to show scientific notation for large numbers
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

# Create legend
legend_elements = [
    mpatches.Patch(color='black', label='Without treatment'),
    mpatches.Patch(color='green', label='Chemotherapy'),
    mpatches.Patch(color='red', label='Radiotherapy'),
    mpatches.Patch(color='magenta', label='Immunotherapy')
]
fig.legend(handles=legend_elements, loc='lower center', ncol=4, fontsize=11, 
          bbox_to_anchor=(0.5, -0.02), frameon=True)

plt.tight_layout(rect=[0, 0.03, 1, 1])
plt.savefig('data/protocols/Figure_4_complete.png', dpi=300, bbox_inches='tight')
print("Figure saved to: data/protocols/Figure_4_complete.png")
plt.show()