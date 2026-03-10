import os
import shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
from protocols import get_params, get_default_ic, run_simulation, get_param_index

def concatenate_simulations(df_list):
    """Concatenate multiple simulation dataframes, adjusting time values"""
    if not df_list:
        return None
    
    combined_df = df_list[0].copy()
    current_time = df_list[0]['time'].iloc[-1]
    
    for df in df_list[1:]:
        df_shifted = df.copy()
        df_shifted['time'] = df_shifted['time'] + current_time
        combined_df = pd.concat([combined_df, df_shifted.iloc[1:]], ignore_index=True)
        current_time = combined_df['time'].iloc[-1]
    
    return combined_df

def add_phase_colors(df, phase_dict):
    """Add phase information and colors to dataframe"""
    colors = []
    for t in df['time']:
        color = 'black'  # default
        for (t_start, t_end), phase_color in phase_dict.items():
            if t_start <= t < t_end:
                color = phase_color
                break
        colors.append(color)
    df['color'] = colors
    return df

def protocol_1_without_resistant(mC=0, mS=0):
    """Protocol 1: Without Resistant Cells"""
    params = get_params()
    params[get_param_index('mC')] = mC
    params[get_param_index('mS')] = mS
    ic = get_default_ic()
    
    print("Protocol 1 (No Resistant Cells) - Running simulations...")
    
    print("  Phase 1: Growth (0-200 days)")
    df_growth = run_simulation(ic, params, t_final=200, n_points=500)
    ic_at_200 = df_growth[['S', 'SR', 'C', 'CR', 'M1', 'M2', 'TH1', 'TH2', 'TC', 'Treg', 'IL10', 'IFNgamma', 'IL2']].iloc[-1].values
    
    print("  Phase 2: Chemotherapy (200-280 days)")
    f_c = 0.3  # Reduced chemotherapy frequency
    M_c = 0.03  # Reduced efficiency
    d_c = 10.0
    u2_C_chemo = f_c * (1 - np.exp(-M_c * d_c))
    df_chemo = run_simulation(ic_at_200, params, t_final=80, n_points=400, u2_C=u2_C_chemo)
    ic_after_chemo = df_chemo[['S', 'SR', 'C', 'CR', 'M1', 'M2', 'TH1', 'TH2', 'TC', 'Treg', 'IL10', 'IFNgamma', 'IL2']].iloc[-1].values
    
    print("  Phase 3: Radiotherapy (280-308 days)")
    alpha = 0.15  # Moderate radiation sensitivity
    beta = 0.01   # Moderate dose response
    d_R = 2.5
    u1_radio = 1 - np.exp(-alpha * d_R - beta * d_R**2)
    df_radio = run_simulation(ic_after_chemo, params, t_final=28, n_points=200, u1=u1_radio)
    ic_after_radio = df_radio[['S', 'SR', 'C', 'CR', 'M1', 'M2', 'TH1', 'TH2', 'TC', 'Treg', 'IL10', 'IFNgamma', 'IL2']].iloc[-1].values
    
    print("  Phase 4: Treatment-free (308-500 days)")
    df_free = run_simulation(ic_after_radio, params, t_final=192, n_points=400)
    
    result_df = concatenate_simulations([df_growth, df_chemo, df_radio, df_free])
    
    # Add phase colors
    phase_dict = {
        (0, 200): 'black',
        (200, 280): 'green',
        (280, 308): 'red',
        (308, 500): 'black'
    }
    result_df = add_phase_colors(result_df, phase_dict)
    
    return result_df

def protocol_1_with_resistant(mC=0.01, mS=4e-7):
    """Protocol 1: With Resistant Cells"""
    params = get_params()
    params[get_param_index('mC')] = mC
    params[get_param_index('mS')] = mS
    ic = get_default_ic()
    
    print("Protocol 1 (With Resistant Cells) - Running simulations...")
    
    print("  Phase 1: Growth with resistance (0-200 days)")
    df_growth = run_simulation(ic, params, t_final=200, n_points=500)
    ic_at_200 = df_growth[['S', 'SR', 'C', 'CR', 'M1', 'M2', 'TH1', 'TH2', 'TC', 'Treg', 'IL10', 'IFNgamma', 'IL2']].iloc[-1].values
    
    print("  Phase 2: Chemotherapy (200-280 days)")
    f_c = 0.3
    M_c = 0.03
    d_c = 10.0
    u2_C_chemo = f_c * (1 - np.exp(-M_c * d_c))
    df_chemo = run_simulation(ic_at_200, params, t_final=80, n_points=400, u2_C=u2_C_chemo)
    ic_after_chemo = df_chemo[['S', 'SR', 'C', 'CR', 'M1', 'M2', 'TH1', 'TH2', 'TC', 'Treg', 'IL10', 'IFNgamma', 'IL2']].iloc[-1].values
    
    print("  Phase 3: Radiotherapy (280-308 days)")
    alpha = 0.15
    beta = 0.01
    d_R = 2.5
    u1_radio = 1 - np.exp(-alpha * d_R - beta * d_R**2)
    df_radio = run_simulation(ic_after_chemo, params, t_final=28, n_points=200, u1=u1_radio)
    ic_after_radio = df_radio[['S', 'SR', 'C', 'CR', 'M1', 'M2', 'TH1', 'TH2', 'TC', 'Treg', 'IL10', 'IFNgamma', 'IL2']].iloc[-1].values
    
    print("  Phase 4: Treatment-free (308-500 days)")
    df_free = run_simulation(ic_after_radio, params, t_final=192, n_points=400)
    
    result_df = concatenate_simulations([df_growth, df_chemo, df_radio, df_free])
    
    # Add phase colors
    phase_dict = {
        (0, 200): 'black',
        (200, 280): 'green',
        (280, 308): 'red',
        (308, 500): 'magenta'
    }
    result_df = add_phase_colors(result_df, phase_dict)
    
    return result_df

def protocol_2_with_resistant(mC=0.01, mS=4e-7):
    """Protocol 2: Combined Therapy"""
    params = get_params()
    params[get_param_index('mC')] = mC
    params[get_param_index('mS')] = mS
    ic = get_default_ic()
    
    print("Protocol 2 (Combined Therapy) - Running simulations...")
    
    print("  Phase 1: Growth with resistance (0-200 days)")
    df_growth = run_simulation(ic, params, t_final=200, n_points=500)
    ic_at_200 = df_growth[['S', 'SR', 'C', 'CR', 'M1', 'M2', 'TH1', 'TH2', 'TC', 'Treg', 'IL10', 'IFNgamma', 'IL2']].iloc[-1].values
    
    print("  Phase 2: Combined therapy (200-240 days)")
    f_c = 0.3
    M_c = 0.03
    d_c = 10.0
    u2_C_combined = f_c * (1 - np.exp(-M_c * d_c))
    
    alpha = 0.15
    beta = 0.01
    d_R = 2.5
    u1_combined = 1 - np.exp(-alpha * d_R - beta * d_R**2)
    
    d_I = 0.3
    M_Tc = 300
    M_TH1 = 300
    u3_Tc_combined = d_I * M_Tc
    u3_TH1_combined = d_I * M_TH1
    
    df_combined = run_simulation(ic_at_200, params, t_final=40, n_points=300,
                                u1=u1_combined, u2_C=u2_C_combined,
                                u3_Tc=u3_Tc_combined, u3_TH1=u3_TH1_combined)
    ic_after_combined = df_combined[['S', 'SR', 'C', 'CR', 'M1', 'M2', 'TH1', 'TH2', 'TC', 'Treg', 'IL10', 'IFNgamma', 'IL2']].iloc[-1].values
    
    print("  Phase 3: Radiotherapy (240-268 days)")
    df_radio = run_simulation(ic_after_combined, params, t_final=28, n_points=200, u1=u1_combined)
    ic_after_radio = df_radio[['S', 'SR', 'C', 'CR', 'M1', 'M2', 'TH1', 'TH2', 'TC', 'Treg', 'IL10', 'IFNgamma', 'IL2']].iloc[-1].values
    
    print("  Phase 4: Treatment-free (268-278 days)")
    df_free = run_simulation(ic_after_radio, params, t_final=10, n_points=100)
    ic_after_free = df_free[['S', 'SR', 'C', 'CR', 'M1', 'M2', 'TH1', 'TH2', 'TC', 'Treg', 'IL10', 'IFNgamma', 'IL2']].iloc[-1].values
    
    print("  Phase 5: Chemotherapy (278-358 days)")
    df_chemo_2 = run_simulation(ic_after_free, params, t_final=80, n_points=400, u2_C=u2_C_combined)
    ic_after_chemo_2 = df_chemo_2[['S', 'SR', 'C', 'CR', 'M1', 'M2', 'TH1', 'TH2', 'TC', 'Treg', 'IL10', 'IFNgamma', 'IL2']].iloc[-1].values
    
    print("  Phase 6: Immunotherapy (358-378 days)")
    d_I_immuno = 0.4
    u3_Tc_immuno = d_I_immuno * M_Tc
    u3_TH1_immuno = d_I_immuno * M_TH1
    df_immuno = run_simulation(ic_after_chemo_2, params, t_final=20, n_points=150,
                              u3_Tc=u3_Tc_immuno, u3_TH1=u3_TH1_immuno)
    
    result_df = concatenate_simulations([df_growth, df_combined, df_radio, df_free, df_chemo_2, df_immuno])
    
    # Add phase colors
    phase_dict = {
        (0, 200): 'black',
        (200, 240): 'green',
        (240, 268): 'red',
        (268, 278): 'black',
        (278, 358): 'green',
        (358, 378): 'magenta'
    }
    result_df = add_phase_colors(result_df, phase_dict)
    
    return result_df

def plot_protocols(p1_no_resist, p1_resist, p2, save_dir='data/protocols'):
    """Create publication-quality plots matching Figure 4 from the paper"""
    os.makedirs(save_dir, exist_ok=True)
    
    fig, axes = plt.subplots(4, 3, figsize=(15, 14))
    
    protocol_data = [
        (p1_no_resist, 'Protocol 1\n(without resistant cells)', 0),
        (p1_resist, 'Protocol 1\n(with resistant cells)', 1),
        (p2, 'Protocol 2\n(with resistant cells)', 2)
    ]
    
    species_info = [
        ('S', 'Stem Cells', 0),
        ('SR', 'Stem Resistant Cells', 1),
        ('C', 'Cancer Cells', 2),
        ('CR', 'Cancer Resistant Cells', 3)
    ]
    
    for df, protocol_title, col_idx in protocol_data:
        for species_name, species_label, row_idx in species_info:
            ax = axes[row_idx, col_idx]
            
            # Get data
            time = df['time'].values
            y_data = df[species_name].values
            colors_array = df['color'].values
            
            # Plot with color transitions
            for i in range(len(time) - 1):
                t_seg = time[i:i+2]
                y_seg = y_data[i:i+2]
                color = colors_array[i]
                ax.plot(t_seg, y_seg, color=color, linewidth=2.5, solid_capstyle='round')
            
            # Formatting
            ax.set_ylabel(species_label, fontsize=11, fontweight='bold')
            ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
            ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            
            if row_idx == 0:
                ax.set_title(protocol_title, fontsize=12, fontweight='bold', pad=10)
            
            if row_idx == 3:
                ax.set_xlabel('Time (days)', fontsize=10, fontweight='bold')
            
            # Set x-axis limits
            if col_idx == 2:  # Protocol 2 extends longer
                ax.set_xlim(0, 800)
            else:
                ax.set_xlim(0, 500)
            
            ax.set_ylim(bottom=0)
            
            # Clean up ticks
            ax.tick_params(labelsize=9)
    
    # Add shared legend at bottom
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color='black', lw=2.5, label='No treatment'),
        Line2D([0], [0], color='green', lw=2.5, label='Chemotherapy'),
        Line2D([0], [0], color='red', lw=2.5, label='Radiotherapy'),
        Line2D([0], [0], color='magenta', lw=2.5, label='Immunotherapy')
    ]
    fig.legend(handles=legend_elements, loc='lower center', ncol=4, fontsize=11,
              bbox_to_anchor=(0.5, -0.01), frameon=True, fancybox=True, shadow=False)
    
    plt.tight_layout(rect=[0, 0.02, 1, 1])
    plt.savefig(f'{save_dir}/Figure_4_Protocols.png', dpi=300, bbox_inches='tight')
    print(f"\n✓ Figure saved: {save_dir}/Figure_4_Protocols.png")
    plt.close()

if __name__ == '__main__':
    os.makedirs('data/protocols', exist_ok=True)
    
    # Run protocols
    print("\n" + "="*70)
    print("PROTOCOL 1: WITHOUT RESISTANT CELLS")
    print("="*70)
    df_p1_no_resist = protocol_1_without_resistant(mC=0, mS=0)
    df_p1_no_resist.to_csv('data/protocols/protocol_1_without_resistant.csv', index=False)
    print("✓ Saved: data/protocols/protocol_1_without_resistant.csv\n")
    
    print("\n" + "="*70)
    print("PROTOCOL 1: WITH RESISTANT CELLS")
    print("="*70)
    df_p1_resist = protocol_1_with_resistant(mC=0.01, mS=4e-7)
    df_p1_resist.to_csv('data/protocols/protocol_1_with_resistant.csv', index=False)
    print("✓ Saved: data/protocols/protocol_1_with_resistant.csv\n")
    
    print("\n" + "="*70)
    print("PROTOCOL 2: COMBINED THERAPY")
    print("="*70)
    df_p2 = protocol_2_with_resistant(mC=0.01, mS=4e-7)
    df_p2.to_csv('data/protocols/protocol_2_combined.csv', index=False)
    print("✓ Saved: data/protocols/protocol_2_combined.csv\n")
    
    # Create plots
    print("\n" + "="*70)
    print("CREATING PLOTS")
    print("="*70)
    plot_protocols(df_p1_no_resist, df_p1_resist, df_p2)
    
    print("\n" + "="*70)
    print("ALL PROTOCOLS COMPLETED SUCCESSFULLY")
    print("="*70 + "\n")
    
    if os.path.exists('__pycache__'):
        shutil.rmtree('__pycache__')