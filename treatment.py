# treatment_updated.py - DOSE-BASED VERSION
# Treatment protocols using d_R, d_c, d_I as primary control parameters

import os
import shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from protocols import (
    get_params, get_default_ic, run_simulation, get_param_index,
    calculate_u1_radiotherapy, calculate_u2S_chemotherapy_stem,
    calculate_u2C_chemotherapy_cancer, calculate_u3Tc_immunotherapy,
    calculate_u3TH1_immunotherapy
)

# ============================================================================
# PROTOCOL STRUCTURE PARAMETERS (Based on Paper Notation)
# ============================================================================
# General Form: DT₂₀₀ → (R^d_R)_t_R → FT_t_FT → (Ch^d_C)_t_C → (I^d_I)_t_I
# Where: DT=Detection Time, R=Radiotherapy, FT=Free Time, Ch=Chemotherapy, I=Immunotherapy
# Superscripts (d) = dosage, Subscripts (t) = duration, (n) = number of cycles/phases
# ============================================================================

# ===== DETECTION PHASE (Protocol baseline) =====
DETECTION_TIME = 200                # Days before any treatment begins (standard across all protocols)

# ===== RADIOTHERAPY PHASE PARAMETERS =====
# Radiotherapy formula: u₁ = 1 - e^(-α·d_R - β·d_R²) where α=0.3, β=0.03
# Each dose d_R (in Gray) kills cancer cells according to linear-quadratic model
GLOBAL_d_R = 60/28                  # **DEFAULT:** 60 Gy / 28 fractions ≈ 2.14 Gy per fraction
                                     # (Paper Protocol 1 standard)
RADIOTHERAPY_PHASE_DURATION = 40    # Duration of single radiotherapy phase: 40 days
NUM_RADIOTHERAPY_PHASES = 1         # Number of radiotherapy phases in Protocol 1 & 2

# ===== TREATMENT-FREE PERIOD =====
FREE_TIME_DURATION = 15             # Days between radiotherapy and second chemotherapy phase
                                     # Allows immune recovery and tumor reoxygenation

# ===== CHEMOTHERAPY PHASE PARAMETERS =====
# Chemotherapy formula: u₂_C = f_c·(1 - e^(-M_C·d_c)) where f_c=0.071, M_C=0.5
# Each dose d_c (in mg/m²) affects cancer cells and stem cells (stem cells have k_S=0.3 protection)
GLOBAL_d_c = 800.0                  # **DEFAULT:** 800 mg/m² per cycle
                                     # (Paper Protocol 1 & 2 standard)
CHEMO_CYCLE_DURATION = 14           # Days per chemotherapy cycle (each dose application)
CHEMO_CYCLES_PER_PHASE = 6          # Number of cycles per chemotherapy phase
CHEMOTHERAPY_PHASE_DURATION = (CHEMO_CYCLE_DURATION * 
                               CHEMO_CYCLES_PER_PHASE)  # = 84 days
NUM_CHEMOTHERAPY_PHASES = 2         # Two chemotherapy phases in Protocol 1 & 2

# ===== IMMUNOTHERAPY PHASE PARAMETERS =====
# Immunotherapy formula: u₃_Tc = d_I·M_Tc, u₃_TH1 = d_I·M_TH1 (M_Tc=0.8, M_TH1=0.7)
# Boosts TC and TH1 cells to enhance anti-tumor immunity
GLOBAL_d_I = 2.0                    # **DEFAULT:** 2.0 immunotherapy dose units
                                     # (Paper Protocol 2 standard)
IMMUNOTHERAPY_PHASE_DURATION = 20   # Duration of immunotherapy phase: 20 days
NUM_IMMUNOTHERAPY_PHASES = 1        # Number of immunotherapy phases (Protocol 2 only)

# ===== OBSERVATION/MONITORING PERIOD =====
# After active treatment ends, system evolves naturally with mutations but no treatment
OBSERVATION_PHASE_START_PROTOCOL1 = 423  # Day treatment ends for Protocol 1
OBSERVATION_PHASE_START_PROTOCOL2 = 443  # Day treatment ends for Protocol 2
FINAL_OBSERVATION_TIME = 800        # Total simulation duration (days)
                                     # Allows observation of post-treatment stabilization/rebound

# ===== TUMOR HETEROGENEITY PARAMETERS =====
GLOBAL_mC = 0.01                    # **DEFAULT:** Cancer cell mutation rate = 1%
                                     # Emergence of drug-resistant cancer cells (CR)
GLOBAL_mS = 4e-7                    # **DEFAULT:** Stem cell mutation rate = 4e-7
                                     # Emergence of drug-resistant stem cells (SR)

# ============================================================================
# CONTROL VARIABLE OVERRIDE FLAG
# ============================================================================
USE_GLOBAL_CONTROL_VARS = True      # If True, calculates u_i from doses
                                     # If False, would use hardcoded values (not recommended)

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def apply_control_variable_overrides(d_R, d_c_S, d_c_C, d_I_Tc, d_I_TH1, phase_name=""):
    """
    Calculate control variables from doses
    
    Parameters:
    -----------
    d_R : float
        Radiotherapy dose (Gray)
    d_c_S, d_c_C : float
        Chemotherapy dose on stem and cancer cells (mg m⁻²)
    d_I_Tc, d_I_TH1 : float
        Immunotherapy doses
    phase_name : str
        Name of treatment phase (for logging)
    
    Returns:
    --------
    tuple
        (u1, u2_S, u2_C, u3_Tc, u3_TH1)
    """
    
    if USE_GLOBAL_CONTROL_VARS:
        # Calculate from doses
        u1 = calculate_u1_radiotherapy(d_R)
        u2_S = calculate_u2S_chemotherapy_stem(d_c_S)
        u2_C = calculate_u2C_chemotherapy_cancer(d_c_C)
        u3_Tc = calculate_u3Tc_immunotherapy(d_I_Tc)
        u3_TH1 = calculate_u3TH1_immunotherapy(d_I_TH1)
        
        print(f"    [{phase_name}] Doses: d_R={d_R:.3f}Gy, d_c={d_c_C:.1f}mg, d_I={d_I_Tc:.2f}")
        print(f"    [{phase_name}] Controls: u1={u1:.4f}, u2_S={u2_S:.4f}, u2_C={u2_C:.4f}, u3_Tc={u3_Tc:.4f}, u3_TH1={u3_TH1:.4f}")
        
        return u1, u2_S, u2_C, u3_Tc, u3_TH1
    else:
        # Direct values (legacy)
        return d_R, d_c_S, d_c_C, d_I_Tc, d_I_TH1

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
# PROTOCOL DEFINITIONS
# ============================================================================

def protocol_0_no_treatment(mC=0.01, mS=4e-7, t_final=800):
    """
    Protocol 0: No Treatment - Pure tumor growth with resistance
    
    Duration: 0-t_final days (typically 800)
    """
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

def protocol_1_chemoradio_free_chemo(mC=0.01, mS=4e-7, t_final=None):
    """
    Protocol 1: Chemo-Radio-Free-Chemo with Post-Treatment Observation
    
    Structure: DT₂₀₀ → (Ch^d_c)_t_c → (R^d_R)_t_R → FT_t_FT → (Ch^d_c)_t_c → OBSERVATION
    
    Timeline (using protocol parameters):
    - 0 to DETECTION_TIME: Detection/Growth (no treatment, allows tumor to grow to detectable size)
    - DETECTION_TIME to DETECTION_TIME+CHEMO: Chemotherapy phase 1 (kills sensitive cells)
    - +CHEMO to +CHEMO+RADIO: Radiotherapy phase (DNA damage to cancer cells)
    - +RADIO to +RADIO+FREE: Treatment-free period (immune recovery, reoxygenation)
    - +FREE to +FREE+CHEMO: Chemotherapy phase 2 (second round of drug)
    - +CHEMO to t_final: OBSERVATION PHASE (no treatment, shows stabilization/rebound)
    
    Key Features:
    - Control variables: u2_S, u2_C active during chemotherapy
    - Control variables: u1 active during radiotherapy
    - Durations parameterized for easy variation in steady-state studies
    """
    if t_final is None:
        t_final = FINAL_OBSERVATION_TIME
    
    params = get_params()
    params[get_param_index('mC')] = mC
    params[get_param_index('mS')] = mS
    ic = get_default_ic()
    
    print(f"Protocol 1 (Chemo-Radio-Free-Chemo) - Running simulations...")
    print(f"  Structure: DT({DETECTION_TIME}d) → Ch({CHEMOTHERAPY_PHASE_DURATION}d) → " +
          f"R({RADIOTHERAPY_PHASE_DURATION}d) → FT({FREE_TIME_DURATION}d) → " +
          f"Ch({CHEMOTHERAPY_PHASE_DURATION}d) → OBS({t_final - OBSERVATION_PHASE_START_PROTOCOL1}d)")
    
    time_point = 0
    
    # Phase 1: Detection (0 to DETECTION_TIME days)
    print(f"  Phase 1: Detection/Growth (0-{DETECTION_TIME} days)")
    df_phase1 = run_simulation(ic, params, t_final=DETECTION_TIME, n_points=400)
    ic_at_detection = df_phase1[['S', 'SR', 'C', 'CR', 'M1', 'M2', 'TH1', 'TH2', 'TC', 'Treg', 'IL10', 'IFNgamma', 'IL2']].iloc[-1].values
    time_point += DETECTION_TIME
    
    # Phase 2: Chemotherapy (DETECTION_TIME to DETECTION_TIME + CHEMOTHERAPY_PHASE_DURATION days)
    print(f"  Phase 2: Chemotherapy ({time_point}-{time_point + CHEMOTHERAPY_PHASE_DURATION} days)")
    u1_c, u2_S_c, u2_C_c, u3_Tc_c, u3_TH1_c = apply_control_variable_overrides(
        0, GLOBAL_d_c, GLOBAL_d_c, 0, 0, "Chemo-1"
    )
    df_chemo1 = run_simulation(ic_at_detection, params, t_final=CHEMOTHERAPY_PHASE_DURATION, n_points=300, 
                              u1=u1_c, u2_S=u2_S_c, u2_C=u2_C_c, u3_Tc=u3_Tc_c, u3_TH1=u3_TH1_c)
    ic_at_chemo1 = df_chemo1[['S', 'SR', 'C', 'CR', 'M1', 'M2', 'TH1', 'TH2', 'TC', 'Treg', 'IL10', 'IFNgamma', 'IL2']].iloc[-1].values
    time_point += CHEMOTHERAPY_PHASE_DURATION
    
    # Phase 3: Radiotherapy (time_point to time_point + RADIOTHERAPY_PHASE_DURATION days)
    print(f"  Phase 3: Radiotherapy ({time_point}-{time_point + RADIOTHERAPY_PHASE_DURATION} days)")
    u1_r, u2_S_r, u2_C_r, u3_Tc_r, u3_TH1_r = apply_control_variable_overrides(
        GLOBAL_d_R, 0, 0, 0, 0, "Radio"
    )
    df_radio = run_simulation(ic_at_chemo1, params, t_final=RADIOTHERAPY_PHASE_DURATION, n_points=200, 
                             u1=u1_r, u2_S=u2_S_r, u2_C=u2_C_r, u3_Tc=u3_Tc_r, u3_TH1=u3_TH1_r)
    ic_at_radio = df_radio[['S', 'SR', 'C', 'CR', 'M1', 'M2', 'TH1', 'TH2', 'TC', 'Treg', 'IL10', 'IFNgamma', 'IL2']].iloc[-1].values
    time_point += RADIOTHERAPY_PHASE_DURATION
    
    # Phase 4: Treatment-free (time_point to time_point + FREE_TIME_DURATION days)
    print(f"  Phase 4: Treatment-free ({time_point}-{time_point + FREE_TIME_DURATION} days)")
    df_free = run_simulation(ic_at_radio, params, t_final=FREE_TIME_DURATION, n_points=100)
    ic_at_free = df_free[['S', 'SR', 'C', 'CR', 'M1', 'M2', 'TH1', 'TH2', 'TC', 'Treg', 'IL10', 'IFNgamma', 'IL2']].iloc[-1].values
    time_point += FREE_TIME_DURATION
    
    # Phase 5: Chemotherapy again (time_point to time_point + CHEMOTHERAPY_PHASE_DURATION days)
    print(f"  Phase 5: Chemotherapy ({time_point}-{time_point + CHEMOTHERAPY_PHASE_DURATION} days)")
    u1_c2, u2_S_c2, u2_C_c2, u3_Tc_c2, u3_TH1_c2 = apply_control_variable_overrides(
        0, GLOBAL_d_c, GLOBAL_d_c, 0, 0, "Chemo-2"
    )
    df_chemo2 = run_simulation(ic_at_free, params, t_final=CHEMOTHERAPY_PHASE_DURATION, n_points=300, 
                              u1=u1_c2, u2_S=u2_S_c2, u2_C=u2_C_c2, u3_Tc=u3_Tc_c2, u3_TH1=u3_TH1_c2)
    ic_at_chemo2 = df_chemo2[['S', 'SR', 'C', 'CR', 'M1', 'M2', 'TH1', 'TH2', 'TC', 'Treg', 'IL10', 'IFNgamma', 'IL2']].iloc[-1].values
    time_point += CHEMOTHERAPY_PHASE_DURATION
    
    # Phase 6: Observation phase (time_point to t_final days)
    df_phases = [df_phase1, df_chemo1, df_radio, df_free, df_chemo2]
    if t_final > time_point:
        remaining_days = t_final - time_point
        print(f"  Phase 6: OBSERVATION ({time_point}-{t_final} days, {remaining_days} days post-treatment)")
        print(f"           ⚠ System evolves naturally with mutations (mC={GLOBAL_mC}, mS={GLOBAL_mS}), NO treatment")
        df_observation = run_simulation(ic_at_chemo2, params, t_final=remaining_days, n_points=300)
        df_phases.append(df_observation)
    
    result_df = concatenate_simulations(df_phases)
    
    # Add phase colors
    phase_dict = {
        (0, DETECTION_TIME): 'black',                                                              # Detection
        (DETECTION_TIME, DETECTION_TIME + CHEMOTHERAPY_PHASE_DURATION): 'green',                 # Chemo 1
        (DETECTION_TIME + CHEMOTHERAPY_PHASE_DURATION, 
         DETECTION_TIME + CHEMOTHERAPY_PHASE_DURATION + RADIOTHERAPY_PHASE_DURATION): 'red',      # Radio
        (DETECTION_TIME + CHEMOTHERAPY_PHASE_DURATION + RADIOTHERAPY_PHASE_DURATION,
         DETECTION_TIME + CHEMOTHERAPY_PHASE_DURATION + RADIOTHERAPY_PHASE_DURATION + FREE_TIME_DURATION): 'black',  # Free
        (DETECTION_TIME + CHEMOTHERAPY_PHASE_DURATION + RADIOTHERAPY_PHASE_DURATION + FREE_TIME_DURATION,
         DETECTION_TIME + 2*CHEMOTHERAPY_PHASE_DURATION + RADIOTHERAPY_PHASE_DURATION + FREE_TIME_DURATION): 'green',  # Chemo 2
        (OBSERVATION_PHASE_START_PROTOCOL1, t_final): 'black'  # Observation
    }
    result_df = add_phase_colors(result_df, phase_dict)
    
    return result_df

def protocol_2_chemoradio_immuno(mC=0.01, mS=4e-7, t_final=None):
    """
    Protocol 2: Chemo-Radio-Free-Chemo-Immuno with Post-Treatment Observation
    
    Structure: DT₂₀₀ → (Ch^d_c)_t_c → (R^d_R)_t_R → FT_t_FT → (Ch^d_c)_t_c → (I^d_I)_t_I → OBSERVATION
    
    Timeline (using protocol parameters):
    - 0 to DETECTION_TIME: Detection/Growth
    - +DETECTION to +CHEMO: Chemotherapy phase 1
    - +CHEMO to +RADIO: Radiotherapy 
    - +RADIO to +FREE: Treatment-free period
    - +FREE to +CHEMO: Chemotherapy phase 2
    - +CHEMO to +IMMUNO: Immunotherapy phase (boosts TC and TH1 cells)
    - +IMMUNO to t_final: OBSERVATION PHASE (post-immunotherapy monitoring)
    
    Key Difference from Protocol 1:
    - Adds IMMUNOTHERAPY_PHASE_DURATION to enhance immune control
    - Control variables: u3_Tc, u3_TH1 active during immunotherapy
    - Later observation period to assess immune-mediated control
    """
    if t_final is None:
        t_final = FINAL_OBSERVATION_TIME
    
    params = get_params()
    params[get_param_index('mC')] = mC
    params[get_param_index('mS')] = mS
    ic = get_default_ic()
    
    print(f"Protocol 2 (Chemo-Radio-Free-Chemo-Immuno) - Running simulations...")
    print(f"  Structure: DT({DETECTION_TIME}d) → Ch({CHEMOTHERAPY_PHASE_DURATION}d) → " +
          f"R({RADIOTHERAPY_PHASE_DURATION}d) → FT({FREE_TIME_DURATION}d) → " +
          f"Ch({CHEMOTHERAPY_PHASE_DURATION}d) → I({IMMUNOTHERAPY_PHASE_DURATION}d) → " +
          f"OBS({t_final - OBSERVATION_PHASE_START_PROTOCOL2}d)")
    
    time_point = 0
    
    # Phase 1: Detection (0 to DETECTION_TIME days)
    print(f"  Phase 1: Detection/Growth (0-{DETECTION_TIME} days)")
    df_phase1 = run_simulation(ic, params, t_final=DETECTION_TIME, n_points=400)
    ic_at_detection = df_phase1[['S', 'SR', 'C', 'CR', 'M1', 'M2', 'TH1', 'TH2', 'TC', 'Treg', 'IL10', 'IFNgamma', 'IL2']].iloc[-1].values
    time_point += DETECTION_TIME
    
    # Phase 2: Chemotherapy
    print(f"  Phase 2: Chemotherapy ({time_point}-{time_point + CHEMOTHERAPY_PHASE_DURATION} days)")
    u1_c, u2_S_c, u2_C_c, u3_Tc_c, u3_TH1_c = apply_control_variable_overrides(
        0, GLOBAL_d_c, GLOBAL_d_c, 0, 0, "Chemo-1"
    )
    df_chemo1 = run_simulation(ic_at_detection, params, t_final=CHEMOTHERAPY_PHASE_DURATION, n_points=300, 
                              u1=u1_c, u2_S=u2_S_c, u2_C=u2_C_c, u3_Tc=u3_Tc_c, u3_TH1=u3_TH1_c)
    ic_at_chemo1 = df_chemo1[['S', 'SR', 'C', 'CR', 'M1', 'M2', 'TH1', 'TH2', 'TC', 'Treg', 'IL10', 'IFNgamma', 'IL2']].iloc[-1].values
    time_point += CHEMOTHERAPY_PHASE_DURATION
    
    # Phase 3: Radiotherapy
    print(f"  Phase 3: Radiotherapy ({time_point}-{time_point + RADIOTHERAPY_PHASE_DURATION} days)")
    u1_r, u2_S_r, u2_C_r, u3_Tc_r, u3_TH1_r = apply_control_variable_overrides(
        GLOBAL_d_R, 0, 0, 0, 0, "Radio"
    )
    df_radio = run_simulation(ic_at_chemo1, params, t_final=RADIOTHERAPY_PHASE_DURATION, n_points=200, 
                             u1=u1_r, u2_S=u2_S_r, u2_C=u2_C_r, u3_Tc=u3_Tc_r, u3_TH1=u3_TH1_r)
    ic_at_radio = df_radio[['S', 'SR', 'C', 'CR', 'M1', 'M2', 'TH1', 'TH2', 'TC', 'Treg', 'IL10', 'IFNgamma', 'IL2']].iloc[-1].values
    time_point += RADIOTHERAPY_PHASE_DURATION
    
    # Phase 4: Treatment-free
    print(f"  Phase 4: Treatment-free ({time_point}-{time_point + FREE_TIME_DURATION} days)")
    df_free = run_simulation(ic_at_radio, params, t_final=FREE_TIME_DURATION, n_points=100)
    ic_at_free = df_free[['S', 'SR', 'C', 'CR', 'M1', 'M2', 'TH1', 'TH2', 'TC', 'Treg', 'IL10', 'IFNgamma', 'IL2']].iloc[-1].values
    time_point += FREE_TIME_DURATION
    
    # Phase 5: Chemotherapy again
    print(f"  Phase 5: Chemotherapy ({time_point}-{time_point + CHEMOTHERAPY_PHASE_DURATION} days)")
    u1_c2, u2_S_c2, u2_C_c2, u3_Tc_c2, u3_TH1_c2 = apply_control_variable_overrides(
        0, GLOBAL_d_c, GLOBAL_d_c, 0, 0, "Chemo-2"
    )
    df_chemo2 = run_simulation(ic_at_free, params, t_final=CHEMOTHERAPY_PHASE_DURATION, n_points=300, 
                              u1=u1_c2, u2_S=u2_S_c2, u2_C=u2_C_c2, u3_Tc=u3_Tc_c2, u3_TH1=u3_TH1_c2)
    ic_at_chemo2 = df_chemo2[['S', 'SR', 'C', 'CR', 'M1', 'M2', 'TH1', 'TH2', 'TC', 'Treg', 'IL10', 'IFNgamma', 'IL2']].iloc[-1].values
    time_point += CHEMOTHERAPY_PHASE_DURATION
    
    # Phase 6: Immunotherapy (new in Protocol 2)
    print(f"  Phase 6: Immunotherapy ({time_point}-{time_point + IMMUNOTHERAPY_PHASE_DURATION} days)")
    print(f"           ⚠ Boosting TC and TH1 cells with d_I={GLOBAL_d_I}")
    u1_i, u2_S_i, u2_C_i, u3_Tc_i, u3_TH1_i = apply_control_variable_overrides(
        0, 0, 0, GLOBAL_d_I, GLOBAL_d_I, "Immuno"
    )
    df_immuno = run_simulation(ic_at_chemo2, params, t_final=IMMUNOTHERAPY_PHASE_DURATION, n_points=150, 
                              u1=u1_i, u2_S=u2_S_i, u2_C=u2_C_i, u3_Tc=u3_Tc_i, u3_TH1=u3_TH1_i)
    ic_at_immuno = df_immuno[['S', 'SR', 'C', 'CR', 'M1', 'M2', 'TH1', 'TH2', 'TC', 'Treg', 'IL10', 'IFNgamma', 'IL2']].iloc[-1].values
    time_point += IMMUNOTHERAPY_PHASE_DURATION
    
    # Phase 7: Observation phase (Protocol 2 specific - post-immunotherapy)
    df_phases = [df_phase1, df_chemo1, df_radio, df_free, df_chemo2, df_immuno]
    if t_final > time_point:
        remaining_days = t_final - time_point
        print(f"  Phase 7: OBSERVATION ({time_point}-{t_final} days, {remaining_days} days post-immunotherapy)")
        print(f"           ⚠ System evolves naturally with mutations (mC={GLOBAL_mC}, mS={GLOBAL_mS}), NO treatment")
        df_observation = run_simulation(ic_at_immuno, params, t_final=remaining_days, n_points=300)
        df_phases.append(df_observation)
    
    result_df = concatenate_simulations(df_phases)
    
    # Add phase colors
    phase_dict = {
        (0, DETECTION_TIME): 'black',                                                              # Detection
        (DETECTION_TIME, DETECTION_TIME + CHEMOTHERAPY_PHASE_DURATION): 'green',                 # Chemo 1
        (DETECTION_TIME + CHEMOTHERAPY_PHASE_DURATION, 
         DETECTION_TIME + CHEMOTHERAPY_PHASE_DURATION + RADIOTHERAPY_PHASE_DURATION): 'red',      # Radio
        (DETECTION_TIME + CHEMOTHERAPY_PHASE_DURATION + RADIOTHERAPY_PHASE_DURATION,
         DETECTION_TIME + CHEMOTHERAPY_PHASE_DURATION + RADIOTHERAPY_PHASE_DURATION + FREE_TIME_DURATION): 'black',  # Free
        (DETECTION_TIME + CHEMOTHERAPY_PHASE_DURATION + RADIOTHERAPY_PHASE_DURATION + FREE_TIME_DURATION,
         DETECTION_TIME + 2*CHEMOTHERAPY_PHASE_DURATION + RADIOTHERAPY_PHASE_DURATION + FREE_TIME_DURATION): 'green',  # Chemo 2
        (DETECTION_TIME + 2*CHEMOTHERAPY_PHASE_DURATION + RADIOTHERAPY_PHASE_DURATION + FREE_TIME_DURATION,
         OBSERVATION_PHASE_START_PROTOCOL2): 'magenta',  # Immunotherapy
        (OBSERVATION_PHASE_START_PROTOCOL2, t_final): 'black'  # Observation
    }
    result_df = add_phase_colors(result_df, phase_dict)
    
    return result_df

# ============================================================================
# PLOTTING FUNCTIONS
# ============================================================================

def save_individual_plots(case_dict, save_dir='data/protocols/individual_plots'):
    """Save individual plots for each case and species"""
    os.makedirs(save_dir, exist_ok=True)
    
    species_info = [
        ('S', 'Stem Cells'),
        ('SR', 'Stem Resistant'),
        ('C', 'Cancer Cells'),
        ('CR', 'Cancer Resistant')
    ]
    
    plot_list = []
    
    for case_name, df in case_dict.items():
        for species_name, species_label in species_info:
            fig, ax = plt.subplots(figsize=(12, 6))
            
            time = df['time'].values
            y_data = df[species_name].values
            colors_array = df['color'].values
            
            # Plot with color transitions
            for i in range(len(time) - 1):
                t_seg = time[i:i+2]
                y_seg = y_data[i:i+2]
                color = colors_array[i]
                ax.plot(t_seg, y_seg, color=color, linewidth=2.5, solid_capstyle='round')
            
            ax.set_xlabel('Time (days)', fontsize=12, fontweight='bold')
            ax.set_ylabel(species_label, fontsize=12, fontweight='bold')
            ax.set_title(f'{case_name} - {species_label}', fontsize=13, fontweight='bold')
            ax.grid(True, alpha=0.3)
            ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            ax.set_xlim(0, 800)
            ax.set_ylim(bottom=0)
            
            filename = f"{case_name.replace(':', '').replace(',', '')}__{species_name}.png"
            filepath = f"{save_dir}/{filename}"
            plt.savefig(filepath, dpi=150, bbox_inches='tight')
            print(f"  Saved: {filepath}")
            
            plot_list.append((filename, case_name, species_label, fig, ax))
            plt.close(fig)
    
    return plot_list

def plot_protocols_combined(case_dict, save_dir='data/protocols'):
    """Create combined 4×6 subplot grid"""
    os.makedirs(save_dir, exist_ok=True)
    
    fig, axes = plt.subplots(4, 6, figsize=(24, 16))
    
    case_order = [
        'Case 1: No treatment, No resistant',
        'Case 2: No treatment, With resistant',
        'Case 3: Protocol 1, No resistant',
        'Case 4: Protocol 1, With resistant',
        'Case 5: Protocol 2, No resistant',
        'Case 6: Protocol 2, With resistant'
    ]
    
    species_info = [
        ('S', 'Stem Cells', 0),
        ('SR', 'Resistant Stem', 1),
        ('C', 'Cancer Cells', 2),
        ('CR', 'Resistant Cancer', 3)
    ]
    
    for col_idx, case_name in enumerate(case_order):
        if case_name not in case_dict:
            continue
        df = case_dict[case_name]
        
        for species_name, species_label, row_idx in species_info:
            ax = axes[row_idx, col_idx]
            
            time = df['time'].values
            y_data = df[species_name].values
            colors_array = df['color'].values
            
            for i in range(len(time) - 1):
                t_seg = time[i:i+2]
                y_seg = y_data[i:i+2]
                color = colors_array[i]
                ax.plot(t_seg, y_seg, color=color, linewidth=2.5, solid_capstyle='round')
            
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
    plt.savefig(f'{save_dir}/treatment_protocols_6cases.png', dpi=300, bbox_inches='tight')
    print(f"\n✓ Combined figure saved: {save_dir}/treatment_protocols_6cases.png")
    plt.close()

# ============================================================================
# MAIN EXECUTION
# ============================================================================

if __name__ == '__main__':
    os.makedirs('data/protocols', exist_ok=True)
    
    print("\n" + "="*80)
    print("DOSE-BASED TREATMENT PROTOCOL ANALYSIS")
    print("="*80)
    print(f"\nGlobal Dose Parameters:")
    print(f"  d_R (Radiotherapy): {GLOBAL_d_R:.3f} Gy/fraction")
    print(f"  d_c (Chemotherapy): {GLOBAL_d_c:.1f} mg m⁻²")
    print(f"  d_I (Immunotherapy): {GLOBAL_d_I:.2f} units")
    print(f"\nMutation Rates:")
    print(f"  mC: {GLOBAL_mC}")
    print(f"  mS: {GLOBAL_mS}")
    print("\n" + "="*80 + "\n")
    
    # Run all 6 cases
    case_dict = {}
    
    print("="*80)
    print("CASE 1: NO TREATMENT, NO RESISTANT CELLS")
    print("="*80)
    df_case1 = protocol_0_no_treatment(mC=0, mS=0)
    df_case1.to_csv('data/protocols/case_1_no_treatment_no_resistant.csv', index=False)
    case_dict['Case 1: No treatment, No resistant'] = df_case1
    print("✓ Saved CSV\n")
    
    print("="*80)
    print("CASE 2: NO TREATMENT, WITH RESISTANT CELLS")
    print("="*80)
    df_case2 = protocol_0_no_treatment(mC=GLOBAL_mC, mS=GLOBAL_mS)
    df_case2.to_csv('data/protocols/case_2_no_treatment_resistant.csv', index=False)
    case_dict['Case 2: No treatment, With resistant'] = df_case2
    print("✓ Saved CSV\n")
    
    print("="*80)
    print("CASE 3: PROTOCOL 1, NO RESISTANT CELLS")
    print("="*80)
    df_case3 = protocol_1_chemoradio_free_chemo(mC=0, mS=0)
    df_case3.to_csv('data/protocols/case_3_protocol1_no_resistant.csv', index=False)
    case_dict['Case 3: Protocol 1, No resistant'] = df_case3
    print("✓ Saved CSV\n")
    
    print("="*80)
    print("CASE 4: PROTOCOL 1, WITH RESISTANT CELLS")
    print("="*80)
    df_case4 = protocol_1_chemoradio_free_chemo(mC=GLOBAL_mC, mS=GLOBAL_mS)
    df_case4.to_csv('data/protocols/case_4_protocol1_resistant.csv', index=False)
    case_dict['Case 4: Protocol 1, With resistant'] = df_case4
    print("✓ Saved CSV\n")
    
    print("="*80)
    print("CASE 5: PROTOCOL 2, NO RESISTANT CELLS")
    print("="*80)
    df_case5 = protocol_2_chemoradio_immuno(mC=0, mS=0)
    df_case5.to_csv('data/protocols/case_5_protocol2_no_resistant.csv', index=False)
    case_dict['Case 5: Protocol 2, No resistant'] = df_case5
    print("✓ Saved CSV\n")
    
    print("="*80)
    print("CASE 6: PROTOCOL 2, WITH RESISTANT CELLS")
    print("="*80)
    df_case6 = protocol_2_chemoradio_immuno(mC=GLOBAL_mC, mS=GLOBAL_mS)
    df_case6.to_csv('data/protocols/case_6_protocol2_resistant.csv', index=False)
    case_dict['Case 6: Protocol 2, With resistant'] = df_case6
    print("✓ Saved CSV\n")
    
    # Create combined plot
    print("="*80)
    print("CREATING COMBINED PLOT (4×6 GRID)")
    print("="*80)
    plot_protocols_combined(case_dict)
    
    # Save individual plots
    print("\n" + "="*80)
    print("SAVING INDIVIDUAL PLOTS (24 total)")
    print("="*80)
    plot_list = save_individual_plots(case_dict)
    
    print("\n" + "="*80)
    print("ALL CASES COMPLETED SUCCESSFULLY")
    print("="*80 + "\n")
    
    if os.path.exists('__pycache__'):
        shutil.rmtree('__pycache__')