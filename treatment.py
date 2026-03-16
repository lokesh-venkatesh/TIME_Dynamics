import os
import shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from protocols import get_params, get_default_ic, run_simulation, get_param_index

# ============================================================================
# GLOBAL CONTROL VARIABLES - Updated defaults for new protocols
# ============================================================================
# Set USE_GLOBAL_CONTROL_VARS = True to use global values instead of calculated ones
USE_GLOBAL_CONTROL_VARS = True

# Cell mutation rates
GLOBAL_mC = 0.01             # Cancer cell mutation rate to resistant (default: 0.01 = 1%)
GLOBAL_mS = 4e-7             # Stem cell mutation rate to resistant (default: 4e-7)

# Treatment control variables
GLOBAL_u1 = 0.33            # Radiotherapy: 60 Gy / 28 fractions
GLOBAL_u2_S = 800.0          # Chemotherapy (Stem cells) - 800 mg dose
GLOBAL_u2_C = 800.0          # Chemotherapy (Cancer cells) - 800 mg dose
GLOBAL_u3_Tc = 2.0           # Immunotherapy (TC cells)
GLOBAL_u3_TH1 = 2.0          # Immunotherapy (TH1 cells)

# ============================================================================
# Helper function to apply global overrides
# ============================================================================
def apply_control_variable_overrides(u1, u2_S, u2_C, u3_Tc, u3_TH1, phase_name=""):
    """
    Apply global control variable overrides if enabled.
    Returns the control variables to use (either calculated or global override).
    """
    if USE_GLOBAL_CONTROL_VARS:
        return GLOBAL_u1, GLOBAL_u2_S, GLOBAL_u2_C, GLOBAL_u3_Tc, GLOBAL_u3_TH1
    else:
        return u1, u2_S, u2_C, u3_Tc, u3_TH1

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

# ============================================================================
# PROTOCOL 0: No Treatment (Baseline Growth)
# ============================================================================
def protocol_0_no_treatment(mC=0.01, mS=4e-7, t_final=800):
    """Protocol 0: No Treatment - Pure tumor growth with resistance"""
    params = get_params()
    params[get_param_index('mC')] = mC
    params[get_param_index('mS')] = mS
    ic = get_default_ic()
    
    print(f"Protocol 0 (No Treatment) - Running simulation...")
    print(f"  Growth (0-{t_final} days)")
    df_growth = run_simulation(ic, params, t_final=t_final, n_points=1000)
    
    # Add phase colors (all black - no treatment)
    phase_dict = {(0, t_final): 'black'}
    df_growth = add_phase_colors(df_growth, phase_dict)
    
    return df_growth

# ============================================================================
# PROTOCOL 1: Chemo (84d) → Radio (40d) → Free (15d) → Chemo (84d)
# ============================================================================
def protocol_1_with_resistant(mC=0.01, mS=4e-7, t_final=800):
    """
    Protocol 1: Chemo-Radio-Free-Chemo + No-Treatment Extension
    - 0-200 days: Detection/Growth (black)
    - 200-284 days: Chemotherapy (green) - 14*6 = 84 days
    - 284-324 days: Radiotherapy (red) - 60/28 for 40 days
    - 324-339 days: Treatment-free (black) - 15 days
    - 339-423 days: Chemotherapy (green) - 84 days
    - 423-t_final days: No treatment (black) - for observation
    """
    params = get_params()
    params[get_param_index('mC')] = mC
    params[get_param_index('mS')] = mS
    ic = get_default_ic()
    
    print(f"Protocol 1 (Chemo-Radio-Free-Chemo) - Running simulations...")
    
    # Phase 1: Detection (0-200 days)
    print("  Phase 1: Detection/Growth (0-200 days)")
    df_phase1 = run_simulation(ic, params, t_final=200, n_points=400)
    ic_at_200 = df_phase1[['S', 'SR', 'C', 'CR', 'M1', 'M2', 'TH1', 'TH2', 'TC', 'Treg', 'IL10', 'IFNgamma', 'IL2']].iloc[-1].values
    
    # Phase 2: Chemotherapy (200-284 days, 84 days duration)
    print("  Phase 2: Chemotherapy (200-284 days)")
    u1_c, u2_S_c, u2_C_c, u3_Tc_c, u3_TH1_c = apply_control_variable_overrides(0, 800, 800, 0, 0, "Chemo")
    df_chemo = run_simulation(ic_at_200, params, t_final=84, n_points=300, u1=u1_c, u2_S=u2_S_c, u2_C=u2_C_c, u3_Tc=u3_Tc_c, u3_TH1=u3_TH1_c)
    ic_at_284 = df_chemo[['S', 'SR', 'C', 'CR', 'M1', 'M2', 'TH1', 'TH2', 'TC', 'Treg', 'IL10', 'IFNgamma', 'IL2']].iloc[-1].values
    
    # Phase 3: Radiotherapy (284-324 days, 40 days duration)
    print("  Phase 3: Radiotherapy (284-324 days)")
    u1_r, u2_S_r, u2_C_r, u3_Tc_r, u3_TH1_r = apply_control_variable_overrides(0.355, 0, 0, 0, 0, "Radio")
    df_radio = run_simulation(ic_at_284, params, t_final=40, n_points=200, u1=u1_r, u2_S=u2_S_r, u2_C=u2_C_r, u3_Tc=u3_Tc_r, u3_TH1=u3_TH1_r)
    ic_at_324 = df_radio[['S', 'SR', 'C', 'CR', 'M1', 'M2', 'TH1', 'TH2', 'TC', 'Treg', 'IL10', 'IFNgamma', 'IL2']].iloc[-1].values
    
    # Phase 4: Treatment-free (324-339 days, 15 days duration)
    print("  Phase 4: Treatment-free (324-339 days)")
    df_free = run_simulation(ic_at_324, params, t_final=15, n_points=100)
    ic_at_339 = df_free[['S', 'SR', 'C', 'CR', 'M1', 'M2', 'TH1', 'TH2', 'TC', 'Treg', 'IL10', 'IFNgamma', 'IL2']].iloc[-1].values
    
    # Phase 5: Chemotherapy again (339-423 days, 84 days duration)
    print("  Phase 5: Chemotherapy (339-423 days)")
    u1_c2, u2_S_c2, u2_C_c2, u3_Tc_c2, u3_TH1_c2 = apply_control_variable_overrides(0, 800, 800, 0, 0, "Chemo-2")
    df_chemo2 = run_simulation(ic_at_339, params, t_final=84, n_points=300, u1=u1_c2, u2_S=u2_S_c2, u2_C=u2_C_c2, u3_Tc=u3_Tc_c2, u3_TH1=u3_TH1_c2)
    ic_at_423 = df_chemo2[['S', 'SR', 'C', 'CR', 'M1', 'M2', 'TH1', 'TH2', 'TC', 'Treg', 'IL10', 'IFNgamma', 'IL2']].iloc[-1].values
    
    # Phase 6: No-treatment observation (423 to t_final days)
    df_phases = [df_phase1, df_chemo, df_radio, df_free, df_chemo2]
    if t_final > 423:
        remaining_days = t_final - 423
        print(f"  Phase 6: No-treatment observation ({423}-{t_final} days, {remaining_days} days)")
        df_observation = run_simulation(ic_at_423, params, t_final=remaining_days, n_points=300)
        df_phases.append(df_observation)
    
    result_df = concatenate_simulations(df_phases)
    
    # Add phase colors
    phase_dict = {
        (0, 200): 'black',      # Detection/Growth
        (200, 284): 'green',    # Chemotherapy
        (284, 324): 'red',      # Radiotherapy
        (324, 339): 'black',    # Treatment-free
        (339, 423): 'green',    # Chemotherapy again
        (423, t_final): 'black'  # No-treatment observation
    }
    result_df = add_phase_colors(result_df, phase_dict)
    
    return result_df

# ============================================================================
# PROTOCOL 2: Protocol 1 + Immunotherapy (20d)
# ============================================================================
def protocol_2_with_resistant(mC=0.01, mS=4e-7, t_final=800):
    """
    Protocol 2: Chemo-Radio-Free-Chemo-Immuno + No-Treatment Extension
    Same as Protocol 1, plus:
    - 423-443 days: Immunotherapy (magenta) - 20 days
    - 443-t_final days: No treatment (black) - for observation
    """
    params = get_params()
    params[get_param_index('mC')] = mC
    params[get_param_index('mS')] = mS
    ic = get_default_ic()
    
    print(f"Protocol 2 (Chemo-Radio-Free-Chemo-Immuno) - Running simulations...")
    
    # Phase 1: Detection (0-200 days)
    print("  Phase 1: Detection/Growth (0-200 days)")
    df_phase1 = run_simulation(ic, params, t_final=200, n_points=400)
    ic_at_200 = df_phase1[['S', 'SR', 'C', 'CR', 'M1', 'M2', 'TH1', 'TH2', 'TC', 'Treg', 'IL10', 'IFNgamma', 'IL2']].iloc[-1].values
    
    # Phase 2: Chemotherapy (200-284 days)
    print("  Phase 2: Chemotherapy (200-284 days)")
    u1_c, u2_S_c, u2_C_c, u3_Tc_c, u3_TH1_c = apply_control_variable_overrides(0, 800, 800, 0, 0, "Chemo")
    df_chemo = run_simulation(ic_at_200, params, t_final=84, n_points=300, u1=u1_c, u2_S=u2_S_c, u2_C=u2_C_c, u3_Tc=u3_Tc_c, u3_TH1=u3_TH1_c)
    ic_at_284 = df_chemo[['S', 'SR', 'C', 'CR', 'M1', 'M2', 'TH1', 'TH2', 'TC', 'Treg', 'IL10', 'IFNgamma', 'IL2']].iloc[-1].values
    
    # Phase 3: Radiotherapy (284-324 days)
    print("  Phase 3: Radiotherapy (284-324 days)")
    u1_r, u2_S_r, u2_C_r, u3_Tc_r, u3_TH1_r = apply_control_variable_overrides(0.355, 0, 0, 0, 0, "Radio")
    df_radio = run_simulation(ic_at_284, params, t_final=40, n_points=200, u1=u1_r, u2_S=u2_S_r, u2_C=u2_C_r, u3_Tc=u3_Tc_r, u3_TH1=u3_TH1_r)
    ic_at_324 = df_radio[['S', 'SR', 'C', 'CR', 'M1', 'M2', 'TH1', 'TH2', 'TC', 'Treg', 'IL10', 'IFNgamma', 'IL2']].iloc[-1].values
    
    # Phase 4: Treatment-free (324-339 days)
    print("  Phase 4: Treatment-free (324-339 days)")
    df_free = run_simulation(ic_at_324, params, t_final=15, n_points=100)
    ic_at_339 = df_free[['S', 'SR', 'C', 'CR', 'M1', 'M2', 'TH1', 'TH2', 'TC', 'Treg', 'IL10', 'IFNgamma', 'IL2']].iloc[-1].values
    
    # Phase 5: Chemotherapy again (339-423 days)
    print("  Phase 5: Chemotherapy (339-423 days)")
    u1_c2, u2_S_c2, u2_C_c2, u3_Tc_c2, u3_TH1_c2 = apply_control_variable_overrides(0, 800, 800, 0, 0, "Chemo-2")
    df_chemo2 = run_simulation(ic_at_339, params, t_final=84, n_points=300, u1=u1_c2, u2_S=u2_S_c2, u2_C=u2_C_c2, u3_Tc=u3_Tc_c2, u3_TH1=u3_TH1_c2)
    ic_at_423 = df_chemo2[['S', 'SR', 'C', 'CR', 'M1', 'M2', 'TH1', 'TH2', 'TC', 'Treg', 'IL10', 'IFNgamma', 'IL2']].iloc[-1].values
    
    # Phase 6: Immunotherapy (423-443 days, 20 days duration)
    print("  Phase 6: Immunotherapy (423-443 days)")
    u1_i, u2_S_i, u2_C_i, u3_Tc_i, u3_TH1_i = apply_control_variable_overrides(0, 0, 0, 2.0, 2.0, "Immuno")
    df_immuno = run_simulation(ic_at_423, params, t_final=20, n_points=150, u1=u1_i, u2_S=u2_S_i, u2_C=u2_C_i, u3_Tc=u3_Tc_i, u3_TH1=u3_TH1_i)
    ic_at_443 = df_immuno[['S', 'SR', 'C', 'CR', 'M1', 'M2', 'TH1', 'TH2', 'TC', 'Treg', 'IL10', 'IFNgamma', 'IL2']].iloc[-1].values
    
    # Phase 7: No-treatment observation (443 to t_final days)
    df_phases = [df_phase1, df_chemo, df_radio, df_free, df_chemo2, df_immuno]
    if t_final > 443:
        remaining_days = t_final - 443
        print(f"  Phase 7: No-treatment observation ({443}-{t_final} days, {remaining_days} days)")
        df_observation = run_simulation(ic_at_443, params, t_final=remaining_days, n_points=300)
        df_phases.append(df_observation)
    
    result_df = concatenate_simulations(df_phases)
    
    # Add phase colors
    phase_dict = {
        (0, 200): 'black',      # Detection/Growth
        (200, 284): 'green',    # Chemotherapy
        (284, 324): 'red',      # Radiotherapy
        (324, 339): 'black',    # Treatment-free
        (339, 423): 'green',    # Chemotherapy again
        (423, 443): 'magenta',  # Immunotherapy
        (443, t_final): 'black'  # No-treatment observation
    }
    result_df = add_phase_colors(result_df, phase_dict)
    
    return result_df

# ============================================================================
# Plotting Functions - Save individual plots
# ============================================================================
def save_individual_plots(protocols_dict, save_dir='data/protocols/individual_plots'):
    """
    Save each plot individually (4 species × 5 protocols = 20 plots).
    
    Parameters:
    -----------
    protocols_dict : dict
        Dictionary with protocol names as keys and dataframes as values
    save_dir : str
        Directory to save individual plots
    """
    os.makedirs(save_dir, exist_ok=True)
    
    species_info = [
        ('S', 'Stem Cells', 0),
        ('SR', 'Stem Resistant Cells', 1),
        ('C', 'Cancer Cells', 2),
        ('CR', 'Cancer Resistant Cells', 3)
    ]
    
    plot_list = []
    
    for species_name, species_label, row_idx in species_info:
        for protocol_name, df in protocols_dict.items():
            # Create individual plot
            fig, ax = plt.subplots(figsize=(10, 6))
            
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
            ax.set_xlabel('Time (days)', fontsize=12, fontweight='bold')
            ax.set_ylabel(species_label, fontsize=12, fontweight='bold')
            ax.set_title(f'{protocol_name}: {species_label}', fontsize=14, fontweight='bold', pad=15)
            ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
            ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            ax.tick_params(labelsize=10)
            ax.set_ylim(bottom=0)
            ax.set_xlim(0, 800)
            
            # Add legend
            from matplotlib.lines import Line2D
            legend_elements = [
                Line2D([0], [0], color='black', lw=2.5, label='No treatment'),
                Line2D([0], [0], color='green', lw=2.5, label='Chemotherapy'),
                Line2D([0], [0], color='red', lw=2.5, label='Radiotherapy'),
                Line2D([0], [0], color='magenta', lw=2.5, label='Immunotherapy')
            ]
            ax.legend(handles=legend_elements, loc='upper left', fontsize=10, framealpha=0.95)
            
            plt.tight_layout()
            
            # Save figure
            filename = f'{save_dir}/{protocol_name.replace(" ", "_")}_{species_name}.png'
            plt.savefig(filename, dpi=300, bbox_inches='tight')
            plot_list.append((filename, protocol_name, species_label, fig, ax))
            print(f"✓ Saved: {filename}")
            
    return plot_list

def show_plots_sequentially(plot_list):
    """
    Display plots one by one in a loop (comment out plt.show() if not needed).
    """
    for filename, protocol_name, species_label, fig, ax in plot_list:
        print(f"\nShowing: {protocol_name} - {species_label}")
        # plt.show()  # ← Comment out this line if you don't want to display plots
        plt.close(fig)

def plot_protocols_combined(case_dict, save_dir='data/protocols'):
    """
    Create combined 4×6 subplot grid with 6 cases.
    
    Parameters:
    -----------
    case_dict : dict
        Dictionary with case names and dataframes
    save_dir : str
        Directory to save the figure
    """
    os.makedirs(save_dir, exist_ok=True)
    
    fig, axes = plt.subplots(4, 6, figsize=(24, 16))
    
    # Order of cases for columns
    case_order = [
        'Case 1: No treatment, No resistant',
        'Case 2: No treatment, With resistant',
        'Case 3: Protocol 1, No resistant',
        'Case 4: Protocol 1, With resistant',
        'Case 5: Protocol 2, No resistant',
        'Case 6: Protocol 2, With resistant'
    ]
    
    case_data = [(case_dict[case], case, col_idx) for col_idx, case in enumerate(case_order)]
    
    species_info = [
        ('S', 'Stem Cells', 0),
        ('SR', 'Stem Resistant Cells', 1),
        ('C', 'Cancer Cells', 2),
        ('CR', 'Cancer Resistant Cells', 3)
    ]
    
    for df, case_name, col_idx in case_data:
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
            ax.set_ylabel(species_label, fontsize=10, fontweight='bold')
            ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
            ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            
            if row_idx == 0:
                ax.set_title(case_name, fontsize=11, fontweight='bold', pad=10)
            
            if row_idx == 3:
                ax.set_xlabel('Time (days)', fontsize=10, fontweight='bold')
            
            ax.set_xlim(0, 800)
            ax.set_ylim(bottom=0)
            ax.tick_params(labelsize=8)
    
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
    plt.savefig(f'{save_dir}/Figure_6_Cases_Comparison.png', dpi=300, bbox_inches='tight')
    print(f"\n✓ Combined figure saved: {save_dir}/Figure_6_Cases_Comparison.png")
    plt.close()

if __name__ == '__main__':
    os.makedirs('data/protocols', exist_ok=True)
    
    # Print control variable mode
    print("\n" + "="*80)
    if USE_GLOBAL_CONTROL_VARS:
        print("⚠ USING GLOBAL CONTROL VARIABLE OVERRIDES")
        print(f"  mC={GLOBAL_mC}, mS={GLOBAL_mS}")
        print(f"  u1={GLOBAL_u1}, u2_S={GLOBAL_u2_S}, u2_C={GLOBAL_u2_C}")
        print(f"  u3_Tc={GLOBAL_u3_Tc}, u3_TH1={GLOBAL_u3_TH1}")
    else:
        print("✓ USING CALCULATED CONTROL VARIABLES")
    print("="*80 + "\n")
    
    # Run all 6 cases
    print("="*80)
    print("CASE 1: NO TREATMENT, NO RESISTANT CELLS")
    print("="*80)
    df_case1 = protocol_0_no_treatment(mC=0, mS=0)
    df_case1.to_csv('data/protocols/case_1_no_treatment_no_resistant.csv', index=False)
    print("✓ Saved: data/protocols/case_1_no_treatment_no_resistant.csv\n")
    
    print("="*80)
    print("CASE 2: NO TREATMENT, WITH RESISTANT CELLS")
    print("="*80)
    df_case2 = protocol_0_no_treatment(mC=GLOBAL_mC, mS=GLOBAL_mS)
    df_case2.to_csv('data/protocols/case_2_no_treatment_resistant.csv', index=False)
    print("✓ Saved: data/protocols/case_2_no_treatment_resistant.csv\n")
    
    print("="*80)
    print("CASE 3: PROTOCOL 1, NO RESISTANT CELLS")
    print("="*80)
    df_case3 = protocol_1_with_resistant(mC=0, mS=0)
    df_case3.to_csv('data/protocols/case_3_protocol1_no_resistant.csv', index=False)
    print("✓ Saved: data/protocols/case_3_protocol1_no_resistant.csv\n")
    
    print("="*80)
    print("CASE 4: PROTOCOL 1, WITH RESISTANT CELLS")
    print("="*80)
    df_case4 = protocol_1_with_resistant(mC=GLOBAL_mC, mS=GLOBAL_mS)
    df_case4.to_csv('data/protocols/case_4_protocol1_resistant.csv', index=False)
    print("✓ Saved: data/protocols/case_4_protocol1_resistant.csv\n")
    
    print("="*80)
    print("CASE 5: PROTOCOL 2, NO RESISTANT CELLS")
    print("="*80)
    df_case5 = protocol_2_with_resistant(mC=0, mS=0)
    df_case5.to_csv('data/protocols/case_5_protocol2_no_resistant.csv', index=False)
    print("✓ Saved: data/protocols/case_5_protocol2_no_resistant.csv\n")
    
    print("="*80)
    print("CASE 6: PROTOCOL 2, WITH RESISTANT CELLS")
    print("="*80)
    df_case6 = protocol_2_with_resistant(mC=GLOBAL_mC, mS=GLOBAL_mS)
    df_case6.to_csv('data/protocols/case_6_protocol2_resistant.csv', index=False)
    print("✓ Saved: data/protocols/case_6_protocol2_resistant.csv\n")
    
    # Create combined plot
    print("="*80)
    print("CREATING COMBINED PLOT (4×6 GRID)")
    print("="*80)
    case_dict = {
        'Case 1: No treatment, No resistant': df_case1,
        'Case 2: No treatment, With resistant': df_case2,
        'Case 3: Protocol 1, No resistant': df_case3,
        'Case 4: Protocol 1, With resistant': df_case4,
        'Case 5: Protocol 2, No resistant': df_case5,
        'Case 6: Protocol 2, With resistant': df_case6
    }
    plot_protocols_combined(case_dict)
    
    # Save and show individual plots
    print("\n" + "="*80)
    print("SAVING INDIVIDUAL PLOTS (24 total)")
    print("="*80)
    plot_list = save_individual_plots(case_dict)
    
    # Show plots one by one (comment out plt.show() in show_plots_sequentially() if not needed)
    print("\n" + "="*80)
    print("DISPLAYING INDIVIDUAL PLOTS SEQUENTIALLY")
    print("="*80)
    print("(Uncomment plt.show() in show_plots_sequentially() if you want to see them)\n")
    # show_plots_sequentially(plot_list)  # ← Uncomment to display plots one by one
    
    print("\n" + "="*80)
    print("ALL CASES COMPLETED SUCCESSFULLY")
    print("="*80 + "\n")
    
    if os.path.exists('__pycache__'):
        shutil.rmtree('__pycache__')