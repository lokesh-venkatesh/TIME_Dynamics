# protocols.py - UPDATED VERSION
# Tumor microenvironment ODE model with treatment protocols
# Control variables now calculated from doses (d_R, d_c, d_I)

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

# ============================================================================
# TREATMENT PARAMETER CONSTANTS (From References)
# ============================================================================

# ===== RADIOTHERAPY PARAMETERS (Powathil et al. 2013, Ref 57) =====
ALPHA_RADIATION = 0.3           # Gy⁻¹ - Linear coefficient for LQ model
BETA_RADIATION = 0.03           # Gy⁻² - Quadratic coefficient for LQ model

# Cell-cycle phase sensitivity factors for radiotherapy
SENSITIVITY_SGM = 1.0           # S-G2-M phase (most radiosensitive)
SENSITIVITY_G1 = 0.5            # G1 phase (intermediate)
SENSITIVITY_RESTING = 0.25      # Resting/quiescent phase (least sensitive)

# Oxygen Enhancement Ratio parameters
OER_MAX = 3.0                   # Maximum oxygen enhancement ratio
KM_OXYGEN = 3.0                 # Half-max oxygen concentration (mmHg)

# ===== CHEMOTHERAPY PARAMETERS =====
# f_c: Frequency of chemotherapy application (day⁻¹)
# For protocol: 6 cycles of 14 days each = 84 days total
# Assuming exponential decay model: f_c relates to drug clearance rate
F_C = 0.071                     # ~1/14 days⁻¹ (once per 14-day cycle)

# M_C: Efficiency/sensitivity of chemotherapy (m² mg⁻¹)
# Normalized value - higher means more effective
M_C = 0.5                       # Moderate efficiency (range: 0.1-1.0)

# k_S: Inhibitory effect of IL-4 on stem cell chemotherapy sensitivity
# IL-4 protects stem cells from chemotherapy (reduces drug effectiveness)
# From immunology literature (Refs 44, 59)
K_S = 0.3                       # 30% reduction in chemo effectiveness on stem cells

# ===== IMMUNOTHERAPY PARAMETERS =====
# M_Tc: Sensitivity/responsiveness of TC cells to immunotherapy (mg day⁻¹)
M_TC = 0.8                      # TC sensitivity (range: 0.5-1.0)

# M_TH1: Sensitivity/responsiveness of TH1 cells to immunotherapy (mg day⁻¹)
M_TH1 = 0.7                     # TH1 sensitivity (range: 0.5-1.0)

# ============================================================================
# CONTROL VARIABLE CALCULATION FUNCTIONS
# ============================================================================

def calculate_u1_radiotherapy(d_R):
    """
    Calculate radiotherapy control variable from dose
    
    u_1 = 1 - e^(-alpha*d_R - beta*d_R²)
    
    Parameters:
    -----------
    d_R : float
        Radiotherapy dose in Gray (Gy)
    
    Returns:
    --------
    float
        Probability of cell death due to radiotherapy (0 to 1)
    """
    if d_R <= 0:
        return 0.0
    u1 = 1.0 - np.exp(-ALPHA_RADIATION * d_R - BETA_RADIATION * (d_R ** 2))
    return np.clip(u1, 0.0, 1.0)

def calculate_u2S_chemotherapy_stem(d_c, use_il4_inhibition=True):
    """
    Calculate chemotherapy control variable for stem cells
    
    u_2S = f_c * (1 - e^(-M_C * d_c)) - k_S
    
    Parameters:
    -----------
    d_c : float
        Chemotherapy drug concentration (mg m⁻²)
    use_il4_inhibition : bool
        Whether to apply IL-4 inhibitory effect (default: True)
    
    Returns:
    --------
    float
        Probability of stem cell death due to chemotherapy
    """
    if d_c <= 0:
        return 0.0
    
    # Base chemotherapy effectiveness
    base_effectiveness = F_C * (1.0 - np.exp(-M_C * d_c))
    
    # Apply IL-4 inhibition for stem cells (reduces chemo effectiveness)
    if use_il4_inhibition:
        u2s = base_effectiveness - K_S
    else:
        u2s = base_effectiveness
    
    return np.clip(u2s, 0.0, 1.0)

def calculate_u2C_chemotherapy_cancer(d_c):
    """
    Calculate chemotherapy control variable for cancer cells
    
    u_2C = f_c * (1 - e^(-M_C * d_c))
    
    Parameters:
    -----------
    d_c : float
        Chemotherapy drug concentration (mg m⁻²)
    
    Returns:
    --------
    float
        Probability of cancer cell death due to chemotherapy
    """
    if d_c <= 0:
        return 0.0
    
    u2c = F_C * (1.0 - np.exp(-M_C * d_c))
    return np.clip(u2c, 0.0, 1.0)

def calculate_u3Tc_immunotherapy(d_I):
    """
    Calculate immunotherapy control variable for TC cells
    
    u_3_Tc = d_I * M_Tc
    
    Parameters:
    -----------
    d_I : float
        Immunotherapy dose (units as specified)
    
    Returns:
    --------
    float
        Boost factor for TC cells
    """
    if d_I <= 0:
        return 0.0
    
    u3_tc = d_I * M_TC
    return np.clip(u3_tc, 0.0, 10.0)  # Allow larger values for immune boost

def calculate_u3TH1_immunotherapy(d_I):
    """
    Calculate immunotherapy control variable for TH1 cells
    
    u_3_TH1 = d_I * M_TH1
    
    Parameters:
    -----------
    d_I : float
        Immunotherapy dose (units as specified)
    
    Returns:
    --------
    float
        Boost factor for TH1 cells
    """
    if d_I <= 0:
        return 0.0
    
    u3_th1 = d_I * M_TH1
    return np.clip(u3_th1, 0.0, 10.0)  # Allow larger values for immune boost

# ============================================================================
# ODE MODEL - TUMOR MICROENVIRONMENT (13 equations)
# ============================================================================

def get_param_index(param_name):
    """Get parameter index by name or return integer index"""
    if isinstance(param_name, int):
        return param_name
    
    param_names = {
        # Production rates (beta)
        'betaM2': 0, 'betaC': 1, 'betaCR': 2, 'betaM1': 3, 'betaTc': 4, 'betaTreg': 5,
        # Decay rates (gamma)
        'gammaC': 6, 'gammaCR': 7, 'gammaM1': 8, 'gammaM2': 9, 'gammaTH1': 10,
        'gammaTH2': 11, 'gammaTc': 12, 'gammaTreg': 13, 'gammaIL10': 14,
        # Death rates (delta)
        'deltaC': 15, 'deltaCR': 16, 'deltaM1': 17, 'deltaM2': 18, 'deltaTH1': 19,
        'deltaTH2': 20, 'deltaTc': 21, 'deltaTreg': 22, 'deltaIL10': 23, 'deltaIFN': 24,
        'deltaIL2': 25, 'deltaS': 26, 'deltaSR': 27,
        # Lambda (interaction rates)
        'lambdaM1': 28, 'lambdaM2': 29, 'lambdaTH1': 30, 'lambdaTH2': 31,
        'lambdaTc': 32, 'lambdaTreg': 33, 'lambdaTreg2': 34, 'lambdaIL10': 35,
        # Mu (competition/inhibition)
        'muC1': 36, 'muC2': 37, 'muCR1': 38, 'muCR2': 39, 'muM1': 40,
        'muM2': 41, 'muTH1': 42, 'muTH2': 43, 'muTcCk1': 44,
        # Carrying capacities and other
        'Cmax': 45, 'CRmax': 46,
        # Drug/radiation effects
        'k1': 47, 'k2': 48, 'k3': 49, 'k4': 50, 'k5': 51, 'k6': 52,
        'k7': 53, 'k8': 54, 'k9': 55, 'k10': 56, 'k_s2': 57,
        'ktc1': 58, 'ktc2': 59, 'ktc3': 60, 'ktc4': 61,
        # Mutation and proliferation
        'mC': 62, 'mS': 63, 'p1': 64, 'p2': 65, 'r1': 66, 'r2': 67,
        'tck': 68, 'muM1Ck2': 69, 'muM2Ck1': 70, 'k_base': 71
    }
    
    return param_names.get(param_name, -1)

def get_params():
    """
    Return default parameter vector (71 parameters)
    """
    params = np.array([
        # Production rates (0-5): betaM2, betaC, betaCR, betaM1, betaTc, betaTreg
        0.2, 0.8, 0.001, 0.3, 0.002, 0.001,
        # Decay rates (6-14): gammaC, gammaCR, gammaM1, gammaM2, gammaTH1, gammaTH2, gammaTc, gammaTreg, gammaIL10
        0.05, 0.05, 0.1, 0.08, 0.12, 0.09, 0.15, 0.08, 0.02,
        # Death rates (15-27): deltaC, deltaCR, deltaM1, deltaM2, deltaTH1, deltaTH2, deltaTc, deltaTreg, deltaIL10, deltaIFN, deltaIL2, deltaS, deltaSR
        0.001, 0.0008, 0.05, 0.04, 0.08, 0.06, 0.1, 0.07, 0.01, 0.05, 0.03, 0.001, 0.0005,
        # Lambda (28-35): lambdaM1, lambdaM2, lambdaTH1, lambdaTH2, lambdaTc, lambdaTreg, lambdaTreg2, lambdaIL10
        0.4, 0.3, 0.5, 0.3, 0.6, 0.4, 0.3, 0.2,
        # Mu (36-44): muC1, muC2, muCR1, muCR2, muM1, muM2, muTH1, muTH2, muTcCk1
        0.0001, 0.00005, 0.00008, 0.00004, 0.1, 0.08, 0.12, 0.09, 0.001,
        # Cmax, CRmax (45-46)
        1e8, 1e7,
        # k values (47-61): k1-k10, k_s2, ktc1-ktc4
        0.5, 0.3, 0.4, 0.35, 0.45, 0.38, 0.2, 0.15, 0.25, 0.1, 0.05, 0.1, 0.08, 0.12, 0.09,
        # Mutation rates (62-63): mC, mS
        0.01, 4e-7,
        # p, r, tck values (64-68): p1, p2, r1, r2, tck
        0.2, 0.15, 0.1, 0.08, 0.05,
        # Additional (69-71): muM1Ck2, muM2Ck1, k_base
        0.15, 0.12, 0.02
    ])
    
    return params

def get_default_ic():
    """
    Return default initial conditions (13 species)
    Species: S, SR, C, CR, M1, M2, TH1, TH2, TC, Treg, IL10, IFNgamma, IL2
    """
    return np.array([1, 0, 0, 0, 85000, 15000, 71000, 12000, 56000, 8000, 0.0085, 0.12, 0.0094])

def rhs(t, y, params, u1=0, u2_S=0, u2_C=0, u3_Tc=0, u3_TH1=0):
    """
    Right-hand side of ODE system for tumor microenvironment model
    
    Parameters:
    -----------
    t : float
        Time
    y : array
        Current state [S, SR, C, CR, M1, M2, TH1, TH2, TC, Treg, IL10, IFNgamma, IL2]
    params : array
        Parameter vector (71 parameters)
    u1 : float
        Radiotherapy control variable (0-1)
    u2_S : float
        Chemotherapy on stem cells (0-1)
    u2_C : float
        Chemotherapy on cancer cells (0-1)
    u3_Tc : float
        Immunotherapy on TC cells
    u3_TH1 : float
        Immunotherapy on TH1 cells
    
    Returns:
    --------
    array
        Derivatives dydt
    """
    
    # Unpack state variables
    S, SR, C, CR, M1, M2, TH1, TH2, TC, Treg, IL10, IFNgamma, IL2 = y
    
    # Ensure non-negative populations
    y = np.maximum(y, 0.0)
    y = np.minimum(y, 1e15)  # Prevent numerical overflow
    
    # Extract parameters (using indices)
    betaM2 = params[0]
    betaC = params[1]
    betaCR = params[2]
    betaM1 = params[3]
    betaTc = params[4]
    betaTreg = params[5]
    
    gammaC = params[6]
    gammaCR = params[7]
    gammaM1 = params[8]
    gammaM2 = params[9]
    gammaTH1 = params[10]
    gammaTH2 = params[11]
    gammaTc = params[12]
    gammaTreg = params[13]
    gammaIL10 = params[14]
    
    deltaC = params[15]
    deltaCR = params[16]
    deltaM1 = params[17]
    deltaM2 = params[18]
    deltaTH1 = params[19]
    deltaTH2 = params[20]
    deltaTc = params[21]
    deltaTreg = params[22]
    deltaIL10 = params[23]
    deltaIFN = params[24]
    deltaIL2 = params[25]
    deltaS = params[26]
    deltaSR = params[27]
    
    lambdaM1 = params[28]
    lambdaM2 = params[29]
    lambdaTH1 = params[30]
    lambdaTH2 = params[31]
    lambdaTc = params[32]
    lambdaTreg = params[33]
    lambdaTreg2 = params[34]
    lambdaIL10 = params[35]
    
    muC1 = params[36]
    muC2 = params[37]
    muCR1 = params[38]
    muCR2 = params[39]
    muM1 = params[40]
    muM2 = params[41]
    muTH1 = params[42]
    muTH2 = params[43]
    muTcCk1 = params[44]
    
    Cmax = params[45]
    CRmax = params[46]
    
    k1 = params[47]
    k2 = params[48]
    k3 = params[49]
    k4 = params[50]
    k5 = params[51]
    k6 = params[52]
    k7 = params[53]
    k8 = params[54]
    k9 = params[55]
    k10 = params[56]
    k_s2 = params[57]
    
    ktc1 = params[58]
    ktc2 = params[59]
    ktc3 = params[60]
    ktc4 = params[61]
    
    mC = params[62]
    mS = params[63]
    p1 = params[64]
    p2 = params[65]
    r1 = params[66]
    r2 = params[67]
    tck = params[68]
    muM1Ck2 = params[69]
    muM2Ck1 = params[70]
    
    # Logarithmic terms for tumor growth saturation
    log_eps = 1e-10
    log_C = np.log(max((Cmax + log_eps) / (C + r1 + log_eps), log_eps))
    log_CR = np.log(max((CRmax + log_eps) / (CR + r1 + log_eps), log_eps))
    
    # ===== STEM CELLS (S) =====
    dS = p1 * S * log_C - (mS * S) - deltaS * S - u2_S * S
    
    # ===== RESISTANT STEM CELLS (SR) =====
    dSR = (mS * S) - deltaSR * SR - u2_S * (1 - k_s2) * SR
    
    # ===== CANCER CELLS (C) =====
    dC = (betaC * S) + (p2 * C * log_C) - (mC * C) - (gammaC * C * IL10 / (muC1 + IL10)) - \
         (k1 * C * TC / (muC2 + TC)) - deltaC * C - u2_C * C
    
    # ===== RESISTANT CANCER CELLS (CR) =====
    dCR = (betaCR * SR) + (mC * C) + (p2 * CR * log_CR) - (gammaCR * CR * IL10 / (muCR1 + IL10)) - \
          (k2 * CR * TC / (muCR2 + TC)) - deltaCR * CR - 0.5 * u2_C * CR
    
    # ===== MACROPHAGE M1 =====
    dM1 = (betaM1) - (lambdaM1 * M1 * IFNgamma / (muM1 + IFNgamma)) - \
          (deltaM1 * M1) - (k3 * M1 * C / (muM1Ck2 + C))
    
    # ===== MACROPHAGE M2 =====
    dM2 = (betaM2) + (lambdaM1 * M1 * IFNgamma / (muM1 + IFNgamma)) - (deltaM2 * M2) - \
          (k4 * M2 * C / (muM2Ck1 + C))
    
    # ===== T HELPER 1 CELLS (TH1) =====
    dTH1 = (lambdaTH1 * TH1 * IL2 / (muTH1 + IL2)) - (k5 * TH1 * IL10 / (muTH1 + IL10)) - \
           (gammaTH1 * TH1) - (deltaTH1 * TH1) + u3_TH1 * 0.1 * TH1
    
    # ===== T HELPER 2 CELLS (TH2) =====
    dTH2 = (lambdaTH2 * TH2 * IL10 / (muTH2 + IL10)) - (k6 * TH2 * IFNgamma / (muTH2 + IFNgamma)) - \
           (gammaTH2 * TH2) - (deltaTH2 * TH2)
    
    # ===== CYTOTOXIC T CELLS (TC) =====
    dTC = (lambdaTc * TC * IL2 / (muTcCk1 + IL2)) - (gammaTc * TC) - (deltaTc * TC) + u3_Tc * 0.1 * TC
    
    # ===== REGULATORY T CELLS (Treg) =====
    dTreg = (lambdaTreg * Treg * IL10 / (muTH1 + IL10)) + (lambdaTreg2 * IL10 / (muTH1 + IL10)) - \
            (k7 * Treg) - (gammaTreg * Treg) - (deltaTreg * Treg)
    
    # ===== CYTOKINES: IL-10 =====
    dIL10 = (lambdaIL10 * M2) + (k8 * Treg) - (k9 * IL10 * TC / (muTcCk1 + TC)) - \
            (gammaIL10 * IL10) - (deltaIL10 * IL10)
    
    # ===== CYTOKINES: IFNgamma =====
    dIFNgamma = (k10 * TH1) - (k5 * IFNgamma * TH2 / (muTH2 + IFNgamma)) - (deltaIFN * IFNgamma)
    
    # ===== CYTOKINES: IL-2 =====
    dIL2 = (tck * TC) + (ktc1 * TH1) - (ktc2 * IL2 * Treg / (muTcCk1 + IL2)) - \
           (ktc3 * IL2) - (ktc4 * IL2 * IL10 / (muTcCk1 + IL10)) - (deltaIL2 * IL2)
    
    # Compile derivatives
    derivs = np.array([dS, dSR, dC, dCR, dM1, dM2, dTH1, dTH2, dTC, dTreg, dIL10, dIFNgamma, dIL2])
    
    # Handle NaNs and Infs
    derivs = np.nan_to_num(derivs, nan=0.0, posinf=0.0, neginf=0.0)
    
    return derivs

def run_simulation(ic, params, t_final=800, n_points=1000, u1=0, u2_S=0, u2_C=0, u3_Tc=0, u3_TH1=0):
    """
    Run ODE simulation
    
    Parameters:
    -----------
    ic : array
        Initial conditions
    params : array
        Parameter vector
    t_final : float
        Final time (days)
    n_points : int
        Number of time points
    u1, u2_S, u2_C, u3_Tc, u3_TH1 : float
        Control variables
    
    Returns:
    --------
    DataFrame
        Simulation results with time and all state variables
    """
    
    t_span = (0, t_final)
    t_eval = np.linspace(0, t_final, n_points)
    
    # Solve ODE
    sol = solve_ivp(rhs, t_span, ic, args=(params, u1, u2_S, u2_C, u3_Tc, u3_TH1),
                    method='BDF', t_eval=t_eval, rtol=1e-10, atol=1e-12, dense_output=True)
    
    # Create output DataFrame
    species_names = ['S', 'SR', 'C', 'CR', 'M1', 'M2', 'TH1', 'TH2', 'TC', 'Treg', 'IL10', 'IFNgamma', 'IL2']
    df = pd.DataFrame(sol.y.T, columns=species_names)
    df['time'] = sol.t
    
    return df[['time'] + species_names]