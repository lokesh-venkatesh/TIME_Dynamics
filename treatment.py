import os
import shutil
import numpy as np
import pandas as pd
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

def protocol_1_without_resistant(mC=0, mS=0):
    """
    Protocol 1: Without Resistant Cells (mC=0, mS=0)
    - Growth phase (0-200 days)
    - Chemotherapy cycles
    - Radiotherapy cycles
    - Treatment-free period
    """
    params = get_params()
    params[get_param_index('mC')] = mC
    params[get_param_index('mS')] = mS
    ic = get_default_ic()
    
    print("Protocol 1 (No Resistant Cells) - Running simulations...")
    
    # Phase 1: Growth until detection (0-200 days)
    print("  Phase 1: Growth (0-200 days)")
    df_growth = run_simulation(ic, params, t_final=200, n_points=500)
    ic_at_200 = df_growth[['S', 'SR', 'C', 'CR', 'M1', 'M2', 'TH1', 'TH2', 'TC', 'Treg', 'IL10', 'IFNgamma', 'IL2']].iloc[-1].values
    
    # Phase 2: Chemotherapy cycles (200-300 days, 50 days treatment)
    print("  Phase 2: Chemotherapy (200-300 days)")
    alpha_chemo = 0.05
    beta_chemo = 0.001
    dose_chemo = 5.0
    u2_C_chemo = alpha_chemo * (1 - np.exp(-beta_chemo * dose_chemo))
    df_chemo = run_simulation(ic_at_200, params, t_final=100, n_points=500, u2_C=u2_C_chemo)
    ic_after_chemo = df_chemo[['S', 'SR', 'C', 'CR', 'M1', 'M2', 'TH1', 'TH2', 'TC', 'Treg', 'IL10', 'IFNgamma', 'IL2']].iloc[-1].values
    
    # Phase 3: Radiotherapy cycles (300-400 days, 50 days treatment)
    print("  Phase 3: Radiotherapy (300-400 days)")
    alpha_radio = 0.3
    beta_radio = 0.05
    dose_radio = 3.0
    u1_radio = 1 - np.exp(-alpha_radio * dose_radio - beta_radio * dose_radio**2)
    df_radio = run_simulation(ic_after_chemo, params, t_final=100, n_points=500, u1=u1_radio)
    ic_after_radio = df_radio[['S', 'SR', 'C', 'CR', 'M1', 'M2', 'TH1', 'TH2', 'TC', 'Treg', 'IL10', 'IFNgamma', 'IL2']].iloc[-1].values
    
    # Phase 4: Treatment-free period (400-500 days)
    print("  Phase 4: Treatment-free (400-500 days)")
    df_free = run_simulation(ic_after_radio, params, t_final=100, n_points=500)
    
    # Combine all phases
    result_df = concatenate_simulations([df_growth, df_chemo, df_radio, df_free])
    
    return result_df

def protocol_1_with_resistant(mC=0.01, mS=4e-7):
    """
    Protocol 1: With Resistant Cells (mC=0.01, mS=4e-7)
    - Growth phase with resistance development (0-200 days)
    - Chemotherapy cycles
    - Radiotherapy cycles
    - Treatment-free period
    """
    params = get_params()
    params[get_param_index('mC')] = mC
    params[get_param_index('mS')] = mS
    ic = get_default_ic()
    
    print("Protocol 1 (With Resistant Cells) - Running simulations...")
    
    # Phase 1: Growth until detection (0-200 days)
    print("  Phase 1: Growth with resistance (0-200 days)")
    df_growth = run_simulation(ic, params, t_final=200, n_points=500)
    ic_at_200 = df_growth[['S', 'SR', 'C', 'CR', 'M1', 'M2', 'TH1', 'TH2', 'TC', 'Treg', 'IL10', 'IFNgamma', 'IL2']].iloc[-1].values
    
    # Phase 2: Chemotherapy cycles (200-300 days)
    print("  Phase 2: Chemotherapy (200-300 days)")
    alpha_chemo = 0.05
    beta_chemo = 0.001
    dose_chemo = 5.0
    u2_C_chemo = alpha_chemo * (1 - np.exp(-beta_chemo * dose_chemo))
    df_chemo = run_simulation(ic_at_200, params, t_final=100, n_points=500, u2_C=u2_C_chemo)
    ic_after_chemo = df_chemo[['S', 'SR', 'C', 'CR', 'M1', 'M2', 'TH1', 'TH2', 'TC', 'Treg', 'IL10', 'IFNgamma', 'IL2']].iloc[-1].values
    
    # Phase 3: Radiotherapy cycles (300-400 days)
    print("  Phase 3: Radiotherapy (300-400 days)")
    alpha_radio = 0.3
    beta_radio = 0.05
    dose_radio = 3.0
    u1_radio = 1 - np.exp(-alpha_radio * dose_radio - beta_radio * dose_radio**2)
    df_radio = run_simulation(ic_after_chemo, params, t_final=100, n_points=500, u1=u1_radio)
    ic_after_radio = df_radio[['S', 'SR', 'C', 'CR', 'M1', 'M2', 'TH1', 'TH2', 'TC', 'Treg', 'IL10', 'IFNgamma', 'IL2']].iloc[-1].values
    
    # Phase 4: Treatment-free period (400-500 days)
    print("  Phase 4: Treatment-free (400-500 days)")
    df_free = run_simulation(ic_after_radio, params, t_final=100, n_points=500)
    
    # Combine all phases
    result_df = concatenate_simulations([df_growth, df_chemo, df_radio, df_free])
    
    return result_df

def protocol_2_with_resistant(mC=0.01, mS=4e-7):
    """
    Protocol 2: Combined Chemo + Radio + Immunotherapy (With Resistant Cells)
    - Growth phase with resistance development (0-200 days)
    - Combined therapy cycles with immunotherapy boost
    - Multiple cycles optimized for synergy
    """
    params = get_params()
    params[get_param_index('mC')] = mC
    params[get_param_index('mS')] = mS
    ic = get_default_ic()
    
    print("Protocol 2 (Combined Therapy) - Running simulations...")
    
    # Phase 1: Growth until detection (0-200 days)
    print("  Phase 1: Growth with resistance (0-200 days)")
    df_growth = run_simulation(ic, params, t_final=200, n_points=500)
    ic_at_200 = df_growth[['S', 'SR', 'C', 'CR', 'M1', 'M2', 'TH1', 'TH2', 'TC', 'Treg', 'IL10', 'IFNgamma', 'IL2']].iloc[-1].values
    
    # Phase 2: Combined chemotherapy + radiotherapy + immunotherapy (200-300 days)
    print("  Phase 2: Combined therapy (200-300 days)")
    alpha_chemo = 0.05
    beta_chemo = 0.001
    dose_chemo = 5.0
    u2_C_chemo = alpha_chemo * (1 - np.exp(-beta_chemo * dose_chemo))
    
    alpha_radio = 0.3
    beta_radio = 0.05
    dose_radio = 3.0
    u1_radio = 1 - np.exp(-alpha_radio * dose_radio - beta_radio * dose_radio**2)
    
    # Immunotherapy boost
    d_I = 2.0
    M_Tc = 500
    M_TH1 = 500
    u3_Tc = d_I * M_Tc
    u3_TH1 = d_I * M_TH1
    
    df_combined = run_simulation(ic_at_200, params, t_final=100, n_points=500, 
                                u1=u1_radio, u2_C=u2_C_chemo, u3_Tc=u3_Tc, u3_TH1=u3_TH1)
    ic_after_combined = df_combined[['S', 'SR', 'C', 'CR', 'M1', 'M2', 'TH1', 'TH2', 'TC', 'Treg', 'IL10', 'IFNgamma', 'IL2']].iloc[-1].values
    
    # Phase 3: Immunotherapy maintenance (300-400 days, reduced therapy)
    print("  Phase 3: Immunotherapy maintenance (300-400 days)")
    d_I_maintenance = 1.5
    u3_Tc_maint = d_I_maintenance * M_Tc
    u3_TH1_maint = d_I_maintenance * M_TH1
    
    df_immuno = run_simulation(ic_after_combined, params, t_final=100, n_points=500,
                              u3_Tc=u3_Tc_maint, u3_TH1=u3_TH1_maint)
    ic_after_immuno = df_immuno[['S', 'SR', 'C', 'CR', 'M1', 'M2', 'TH1', 'TH2', 'TC', 'Treg', 'IL10', 'IFNgamma', 'IL2']].iloc[-1].values
    
    # Phase 4: Treatment-free period (400-500 days)
    print("  Phase 4: Treatment-free (400-500 days)")
    df_free = run_simulation(ic_after_immuno, params, t_final=100, n_points=500)
    
    # Phase 5: Extended observation (500-800 days)
    print("  Phase 5: Extended observation (500-800 days)")
    df_extended = run_simulation(df_free[['S', 'SR', 'C', 'CR', 'M1', 'M2', 'TH1', 'TH2', 'TC', 'Treg', 'IL10', 'IFNgamma', 'IL2']].iloc[-1].values, 
                                params, t_final=300, n_points=500)
    
    # Combine all phases
    result_df = concatenate_simulations([df_growth, df_combined, df_immuno, df_free, df_extended])
    
    return result_df

if __name__ == '__main__':
    os.makedirs('data/protocols', exist_ok=True)
    
    # Run Protocol 1 - Without resistant cells
    print("\n" + "="*60)
    print("PROTOCOL 1: WITHOUT RESISTANT CELLS")
    print("="*60)
    df_p1_no_resist = protocol_1_without_resistant(mC=0, mS=0)
    df_p1_no_resist.to_csv('data/protocols/protocol_1_without_resistant.csv', index=False)
    print("Saved: data/protocols/protocol_1_without_resistant.csv\n")
    
    # Run Protocol 1 - With resistant cells
    print("\n" + "="*60)
    print("PROTOCOL 1: WITH RESISTANT CELLS")
    print("="*60)
    df_p1_resist = protocol_1_with_resistant(mC=0.01, mS=4e-7)
    df_p1_resist.to_csv('data/protocols/protocol_1_with_resistant.csv', index=False)
    print("Saved: data/protocols/protocol_1_with_resistant.csv\n")
    
    # Run Protocol 2 - With resistant cells and combined therapy
    print("\n" + "="*60)
    print("PROTOCOL 2: COMBINED THERAPY (CHEMO + RADIO + IMMUNO)")
    print("="*60)
    df_p2 = protocol_2_with_resistant(mC=0.01, mS=4e-7)
    df_p2.to_csv('data/protocols/protocol_2_combined.csv', index=False)
    print("Saved: data/protocols/protocol_2_combined.csv\n")
    
    print("\n" + "="*60)
    print("ALL PROTOCOLS COMPLETED SUCCESSFULLY")
    print("="*60)
    
    if os.path.exists('__pycache__'):
        shutil.rmtree('__pycache__')