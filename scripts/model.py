import numpy as np
from scipy.integrate import solve_ivp
import pandas as pd

log_eps = 0 #1e-12

PARAM_NAMES = { # Map parameter names to indices
    'betaM2': 0, 'betaTc': 1, 'betaTh1CK2': 2, 'betaTh1CK3': 3, 'betaTh2': 4, 'betaTreg': 5,
    'gammaC': 6, 'gammaCR': 7, 'gammaM1': 8, 'gammaM2': 9, 'gammaS': 10, 'gammaTc': 11, 'gammaTh1': 12, 'gammaTh2': 13, 'gammaTreg': 14,
    'deltaC': 15, 'deltaCk1': 16, 'deltaCk2': 17, 'deltaCk3': 18, 'deltaCR': 19, 'deltaM1': 20, 'deltaM2': 21, 'deltaS': 22, 'deltaTc': 23, 'deltaTh1': 24, 'deltaTh2': 25, 'deltaTreg': 26,
    'lambdaM1': 27, 'lambdaM2': 28, 'lambdaTc1': 29, 'lambdaTc2': 30, 'lambdaTc3': 31, 'lambdaTc4': 32, 'lambdaTh1': 33, 'lambdaTh2': 34, 'lambdaTreg2': 35,
    'muC1': 36, 'muC2': 37, 'muS': 38, 'muSR': 39, 'muTcS': 40, 'muTcTreg': 41, 'muTh1Ck1': 42, 'muTh1Ck3': 43, 'muTregCk1': 44,
    'Cmax': 45, 'CRmax': 46,
    'k1': 47, 'k11': 48, 'k2': 49, 'k3': 50, 'k4': 51, 'k5': 52, 'k6': 53, 'k8': 54, 'k9': 55, 'ktc1': 56, 'ktc2': 57, 'ktc3': 58, 'ktc4': 59,
    'mC': 60, 'mS': 61, 'p1': 62, 'p2': 63, 'r1': 64, 'r2': 65, 'tck': 66, 'muM1Ck2': 67, 'muM2Ck1': 68, 'k7': 69, 'k10': 70
}

# Control parameter indices
CONTROL_NAMES = {
    'alpha': 0,          # Radiotherapy coefficient (cell-state dependent)
    'beta': 1,           # Radiotherapy coefficient (cell-state dependent)
    'f_c': 2,            # Frequency of chemotherapy per day
    'M_c': 3,            # Efficiency of chemotherapy drug (m^2 mg^-1)
    'k_S': 4,            # Inhibitory effect of IL-4 on stem cells
    'M_Tc': 5,           # Sensitivity of Tc cells to immunostimulant (mg day^-1)
    'M_TH1': 6           # Sensitivity of TH1 cells to immunostimulant (mg day^-1)
}

def get_param_index(param): # Convert parameter name or index to index
    if isinstance(param, str):
        if param not in PARAM_NAMES:
            raise ValueError(f"Unknown parameter: {param}")
        return PARAM_NAMES[param]
    return param

def get_control_index(control): # Convert control name to index
    if isinstance(control, str):
        if control not in CONTROL_NAMES:
            raise ValueError(f"Unknown control variable: {control}")
        return CONTROL_NAMES[control]
    return control

def get_params(ranges=None):
    # Default values from file (in exact order used in model.rhs)
    defaults = np.array([
        1e-15, 1e-8, 1e-7, 1e-8, 1e-9, 1e-10,  # beta* (6)
        0.1282, 0.1282, 0.7, 0.01, 0.15, 1.0, 2.0, 2.0, 0.3,  # gamma* (9)
        0.8055, 19.757, 6.1212, 8.664339, 5.37e-5, 1.02, 0.05, 2e-7, 5.2939, 2.0, 2.0, 1.0,  # delta* (12)
        1e8, 1e6, 1e5, 5e5, 5e10, 1e5, 1e5, 1e5, 1e7,  # lambda* (9)
        0.75, 0.9, 0.17, 0.18, 1e-10, 1.5e-5, 1e-9, 0.1245, 1e-7,  # mu* (9)
        1e10, 1e10,  # max* (2)
        10.0, 0.001, 10.0, 2.0531, 3.02, 6.7979, 6.9937, 0.01, 0.001,  # k*, ktc* (9)
        1e9, 1e8, 1e9, 1e9, 0.01, 4e-7, 0.2, 0.05, 0.0001, 1e-5, 0.1, 0.01, 0.01, 0.2, 0.01 # some other params (15)
        ])
    params = defaults.copy()
    # if ranges is not None:
    #     for i, (low_mult, high_mult) in ranges.items():
    #        params[i] *= np.random.uniform(low_mult, high_mult)
    return params

def get_default_controls():
    """
    Default control variable values:
    alpha, beta: Radiotherapy coefficients (cell-state dependent)
    f_c: Frequency of chemotherapy per day
    M_c: Efficiency of chemotherapy drug (m^2 mg^-1)
    k_S: Inhibitory effect of IL-4 on stem cells
    M_Tc: Sensitivity of Tc cells to immunostimulant
    M_TH1: Sensitivity of TH1 cells to immunostimulant
    """
    return np.array([
        0.35,      # alpha - radiotherapy coefficient
        0.015,     # beta - radiotherapy coefficient
        0.5,       # f_c - frequency of chemotherapy per day
        0.02,      # M_c - efficiency of chemotherapy drug
        0.1,       # k_S - inhibitory effect of IL-4 on stem cells
        0.1,       # M_Tc - sensitivity of Tc cells
        0.1        # M_TH1 - sensitivity of TH1 cells
    ])

def get_default_ic():
    # Default ICs (WILL STAY THE SAME FOR ALL RUNS IN THIS STUDY)
    return np.array([1, 0, 0, 0, 85000, 15000, 71000, 12000, 56000, 8000, 0.0085, 0.12, 0.0094])

def compute_radiotherapy_control(d_R, alpha, beta):
    """
    Compute u_1 = 1 - exp(-alpha*d_R - beta*d_R^2)
    Probability of cell death due to radiotherapy
    
    Parameters:
    -----------
    d_R : float
        Dose of radiotherapy in Gray units
    alpha : float
        Linear coefficient (cell-state dependent)
    beta : float
        Quadratic coefficient (cell-state dependent)
    
    Returns:
    --------
    u_1 : float
        Probability of cell death (0 to 1)
    """
    exponent = -alpha * d_R - beta * (d_R ** 2)
    u_1 = 1.0 - np.exp(exponent)
    return np.clip(u_1, 0, 1)  # Ensure probability is bounded [0, 1]

def compute_chemotherapy_controls(f_c, M_c, d_c, k_S, IL4):
    """
    Compute chemotherapy control variables
    u_2S = f_c(1 - exp(-M_c*d_c)) - k_S*(IL4/(IL4+threshold))
    u_2C = f_c(1 - exp(-M_c*d_c))
    
    Parameters:
    -----------
    f_c : float
        Frequency of chemotherapy per day
    M_c : float
        Efficiency of chemotherapy drug (m^2 mg^-1)
    d_c : float
        Concentration of the drug (mg m^-2)
    k_S : float
        Inhibitory effect of IL-4 on stem cells
    IL4 : float
        IL-4 concentration (represented by IL10 in the model)
    
    Returns:
    --------
    u_2S, u_2C : tuple of floats
        Death probabilities for stem cells and cancer cells
    """
    base_effect = f_c * (1.0 - np.exp(-M_c * d_c))
    
    # u_2C: applies to drug-sensitive cancer cells
    u_2C = base_effect
    
    # u_2S: applies to drug-sensitive stem cells, reduced by IL-10/IL-4 effect
    # Using a saturation function for IL-10 inhibition
    IL4_threshold = 0.01  # Threshold parameter for IL-10 effect
    inhibition = k_S * (IL4 / (IL4 + IL4_threshold))
    u_2S = base_effect - inhibition
    
    # Ensure probabilities are bounded [0, 1]
    u_2S = np.clip(u_2S, 0, 1)
    u_2C = np.clip(u_2C, 0, 1)
    
    return u_2S, u_2C

def compute_immunotherapy_controls(d_I, M_Tc, M_TH1):
    """
    Compute immunotherapy control variables
    u_3_Tc = d_I * M_Tc
    u_3_TH1 = d_I * M_TH1
    
    Parameters:
    -----------
    d_I : float
        Dose of immunostimulant
    M_Tc : float
        Sensitivity of Tc cells (mg day^-1)
    M_TH1 : float
        Sensitivity of TH1 cells (mg day^-1)
    
    Returns:
    --------
    u_3_Tc, u_3_TH1 : tuple of floats
        Immunotherapy boost factors
    """
    u_3_Tc = d_I * M_Tc
    u_3_TH1 = d_I * M_TH1
    
    return u_3_Tc, u_3_TH1

def rhs(t, y, params, controls=None, d_R=0.0, d_c=0.0, d_I=0.0):
    """
    Right-hand side of the ODE system with control variables
    
    Parameters:
    -----------
    t : float
        Time
    y : ndarray
        State vector [S, SR, C, CR, M1, M2, TH1, TH2, TC, Treg, IL10, IFNgamma, IL2]
    params : ndarray
        Model parameters
    controls : ndarray, optional
        Control parameters [alpha, beta, f_c, M_c, k_S, M_Tc, M_TH1]
    d_R : float, optional
        Radiotherapy dose (Gray units)
    d_c : float, optional
        Chemotherapy drug concentration (mg m^-2)
    d_I : float, optional
        Immunostimulant dose
    
    Returns:
    --------
    dydt : ndarray
        Rate of change of state variables
    """
    # y = np.maximum(y, log_eps)
    S, SR, C, CR, M1, M2, TH1, TH2, TC, Treg, IL10, IFNgamma, IL2 = y
    
    # Default controls if not provided
    if controls is None:
        controls = get_default_controls()
    
    alpha, beta, f_c, M_c, k_S, M_Tc, M_TH1 = controls
    
    # Unpack params (matching order from file)
    betaM2, betaTc, betaTh1CK2, betaTh1CK3, betaTh2, betaTreg, gammaC, gammaCR, gammaM1, gammaM2, gammaS, gammaTc, gammaTh1, gammaTh2, gammaTreg, deltaC, deltaCk1, deltaCk2, deltaCk3, deltaCR, deltaM1, deltaM2, deltaS, deltaTc, deltaTh1, deltaTh2, deltaTreg, lambdaM1, lambdaM2, lambdaTc1, lambdaTc2, lambdaTc3, lambdaTc4, lambdaTh1, lambdaTh2, lambdaTreg2, muC1, muC2, muS, muSR, muTcS, muTcTreg, muTh1Ck1, muTh1Ck3, muTregCk1, Cmax, CRmax, k1, k11, k2, k3, k4, k5, k6, k8, k9, ktc1, ktc2, ktc3, ktc4, mC, mS, p1, p2, r1, r2, tck, muM1Ck2, muM2Ck1, k7, k10 = params
    
    # ========== COMPUTE CONTROL VARIABLES ==========
    # u_1: Radiotherapy effect on cancer cells
    u_1 = compute_radiotherapy_control(d_R, alpha, beta)
    
    # u_2S, u_2C: Chemotherapy effects on stem and cancer cells
    u_2S, u_2C = compute_chemotherapy_controls(f_c, M_c, d_c, k_S, IL10)
    
    # u_3_Tc, u_3_TH1: Immunotherapy boost
    u_3_Tc, u_3_TH1 = compute_immunotherapy_controls(d_I, M_Tc, M_TH1)
    
    # ========== ODE EQUATIONS WITH CONTROL TERMS ==========
    
    # Stem cell S: affected by chemotherapy (u_2S)
    dS = (gammaS*((1-mS)*(1-p1-p2)))*S - (deltaS+(p2*gammaS)+gammaS*(mS*p1/2))*S - ((muS*S*IFNgamma)/(IFNgamma+k1)) - ((tck*S*TC)/(ktc1+TC)) - u_2S*S
    
    # Resistant stem cell SR: affected by chemotherapy (less sensitive)
    dSR = (gammaS*(1-p1-p2) - (deltaS+(p2*gammaS)))*SR + mS*gammaS*(1-p1/2-p2)*S - ((muSR*SR*IFNgamma)/(k2+IFNgamma)) - ((tck*SR*TC)/(ktc2+TC))
    
    # Drug-sensitive cancer cell C: affected by radiotherapy (u_1) and chemotherapy (u_2C)
    dC = gammaC*(1-mC)*np.log((Cmax+log_eps)/(C+r1+log_eps))*C + gammaS*(p1+p2)*S - deltaC*C - mC*gammaC*C + (muC1*C*IL10)/(IL10+k3) - (muC2*C*IFNgamma)/(IFNgamma+k4) - (tck*C*TC)/(ktc3+TC) - u_1*C - u_2C*C
    
    # Drug-resistant cancer cell CR: affected only by radiotherapy (u_1)
    dCR = gammaCR*CR*np.log((CRmax+log_eps)/(CR+r2+log_eps)) + gammaS*SR*(p1+p2) + mC*gammaCR*C - deltaCR*CR + (muC1*CR*IL10)/(IL10+k5) - (muC2*CR*IFNgamma)/(IFNgamma+k6) - (tck*CR*TC)/(ktc4+TC) - u_1*CR
    
    # Macrophage M1
    dM1 = gammaM1*M1*((C+CR)/(M1+lambdaM1)) - deltaM1*M1 + ((muM1Ck2*M1*IFNgamma)/(IFNgamma+k7))
    
    # Macrophage M2
    dM2 = gammaM2*M2*((C+CR)/(M2+lambdaM2)) - deltaM2*M2 + ((muM2Ck1*M2*IL10)/(IL10+k10))
    
    # Cytotoxic T cell TH1: boosted by immunotherapy (u_3_TH1)
    dTH1 = gammaTh1*((TH1*M1)/(lambdaTh1+TH1)) - deltaTh1*TH1 - ((muTh1Ck1*IL10*TH1)/(IL10+k8)) + ((muTh1Ck3*IL2*TH1)/(IL2+k9)) + u_3_TH1*TH1
    
    # Helper T cell TH2
    dTH2 = gammaTh2*((TH2*M2)/(lambdaTh2+TH2)) - deltaTh2*TH2
    
    # Cytotoxic T cell TC: boosted by immunotherapy (u_3_Tc)
    dTC = gammaTc*TC*((C+CR)/(TC+lambdaTc1)) + gammaTc*((TC*TH1)/(TC+lambdaTc4)) - muTcS*TC*((S+SR)/(TC+lambdaTc2)) - deltaTc*TC - muTcTreg*TC*((Treg)/(lambdaTc3+Treg)) + u_3_Tc*TC
    
    # Regulatory T cell Treg
    dTreg = gammaTreg*((Treg*M2)/(Treg+lambdaTreg2)) - deltaTreg*Treg + ((muTregCk1*IL10*Treg)/(Treg+k11))
    
    # Cytokine IL-10 (anti-inflammatory)
    dIL10 = betaM2*M2 - deltaCk1*IL10 + betaTreg*Treg + betaTh2*TH2
    
    # Cytokine IFN-gamma (pro-inflammatory)
    dIFNgamma = betaTh1CK2*TH1 + betaTc*TC - deltaCk2*IFNgamma
    
    # Cytokine IL-2 (T cell proliferation)
    dIL2 = betaTh1CK3*TH1 - deltaCk3*IL2
    
    return np.array([dS, dSR, dC, dCR, dM1, dM2, dTH1, dTH2, dTC, dTreg, dIL10, dIFNgamma, dIL2])

def event_equil(t, y, params, controls=None, d_R=0.0, d_c=0.0, d_I=0.0): 
    return np.linalg.norm(y) * 0.01

event_equil.terminal = True # Stop integration when event is triggered
event_equil.direction = -1 # Only trigger when the function is decreasing (approaching equilibrium)

def run_simulation(initial_conditions, params, controls=None, d_R=0.0, d_c=0.0, d_I=0.0, t_final='equil', n_points=100000):
    """
    Run ODE simulation with control variables
    
    Parameters:
    -----------
    initial_conditions : ndarray
        Initial state vector
    params : ndarray
        Model parameters
    controls : ndarray, optional
        Control parameters [alpha, beta, f_c, M_c, k_S, M_Tc, M_TH1]
    d_R : float, optional
        Radiotherapy dose (Gray units)
    d_c : float, optional
        Chemotherapy drug concentration (mg m^-2)
    d_I : float, optional
        Immunostimulant dose
    t_final : float or 'equil', optional
        Final time or 'equil' to stop at equilibrium
    n_points : int, optional
        Number of output points
    
    Returns:
    --------
    df : pandas.DataFrame
        Solution with columns: time, S, SR, C, CR, M1, M2, TH1, TH2, TC, Treg, IL10, IFNgamma, IL2
    """
    if controls is None:
        controls = get_default_controls()
    
    t_span = [0, 10000.0] if t_final == 'equil' else [0, t_final]
    events = event_equil if t_final == 'equil' else None
    
    sol = solve_ivp(rhs, 
                    t_span, 
                    initial_conditions, 
                    args=(params, controls, d_R, d_c, d_I), 
                    method='BDF', rtol=1e-10, atol=1e-12, 
                    # method='BDF', rtol=1e-6, atol=1e-8, 
                    # method='BDF', rtol=1e-2, atol=1e-4,  # LESS 'STIFF'
                    events=events, 
                    dense_output=True)
    
    if not sol.success:
        raise RuntimeError(f"Solver failed: {sol.message}")
    
    t_eval = np.linspace(0, sol.t[-1], n_points)
    y_eval = sol.sol(t_eval).T
    
    state_names = ['S', 'SR', 'C', 'CR', 'M1', 'M2', 'TH1', 'TH2', 'TC', 'Treg', 'IL10', 'IFNgamma', 'IL2']
    df = pd.DataFrame(y_eval, columns=state_names)
    df.insert(0, 'time', t_eval)
    
    return df