import numpy as np
from scipy.integrate import solve_ivp
import pandas as pd

log_eps = 0

PARAM_NAMES = {
    'betaM2': 0, 'betaTc': 1, 'betaTh1CK2': 2, 'betaTh1CK3': 3, 'betaTh2': 4, 'betaTreg': 5,
    'gammaC': 6, 'gammaCR': 7, 'gammaM1': 8, 'gammaM2': 9, 'gammaS': 10, 'gammaTc': 11, 'gammaTh1': 12, 'gammaTh2': 13, 'gammaTreg': 14,
    'deltaC': 15, 'deltaCk1': 16, 'deltaCk2': 17, 'deltaCk3': 18, 'deltaCR': 19, 'deltaM1': 20, 'deltaM2': 21, 'deltaS': 22, 'deltaTc': 23, 'deltaTh1': 24, 'deltaTh2': 25, 'deltaTreg': 26,
    'lambdaM1': 27, 'lambdaM2': 28, 'lambdaTc1': 29, 'lambdaTc2': 30, 'lambdaTc3': 31, 'lambdaTc4': 32, 'lambdaTh1': 33, 'lambdaTh2': 34, 'lambdaTreg2': 35,
    'muC1': 36, 'muC2': 37, 'muS': 38, 'muSR': 39, 'muTcS': 40, 'muTcTreg': 41, 'muTh1Ck1': 42, 'muTh1Ck3': 43, 'muTregCk1': 44,
    'Cmax': 45, 'CRmax': 46,
    'k1': 47, 'k11': 48, 'k2': 49, 'k3': 50, 'k4': 51, 'k5': 52, 'k6': 53, 'k8': 54, 'k9': 55, 'ktc1': 56, 'ktc2': 57, 'ktc3': 58, 'ktc4': 59,
    'mC': 60, 'mS': 61, 'p1': 62, 'p2': 63, 'r1': 64, 'r2': 65, 'tck': 66, 'muM1Ck2': 67, 'muM2Ck1': 68, 'k7': 69, 'k10': 70
}

def get_param_index(param):
    if isinstance(param, str):
        if param not in PARAM_NAMES:
            raise ValueError(f"Unknown parameter: {param}")
        return PARAM_NAMES[param]
    return param

def get_params(ranges=None):
    defaults = np.array([
        1e-15, 1e-8, 1e-7, 1e-8, 1e-9, 1e-10,
        0.1282, 0.1282, 0.7, 0.01, 0.15, 1.0, 2.0, 2.0, 0.3,
        0.8055, 19.757, 6.1212, 8.664339, 5.37e-5, 1.02, 0.05, 2e-7, 5.2939, 2.0, 2.0, 1.0,
        1e8, 1e6, 1e5, 5e5, 5e10, 1e5, 1e5, 1e5, 1e7,
        0.75, 0.9, 0.17, 0.18, 1e-10, 1.5e-5, 1e-9, 0.1245, 1e-7,
        1e10, 1e10,
        10.0, 0.001, 10.0, 2.0531, 3.02, 6.7979, 6.9937, 0.01, 0.001,
        1e9, 1e8, 1e9, 1e9, 0.01, 4e-7, 0.2, 0.05, 0.0001, 1e-5, 0.1, 0.01, 0.01, 0.2, 0.01
        ])
    params = defaults.copy()
    return params

def get_default_ic():
    return np.array([1, 0, 0, 0, 85000, 15000, 71000, 12000, 56000, 8000, 0.0085, 0.12, 0.0094])

def rhs(t, y, params, u1=0, u2_S=0, u2_C=0, u3_Tc=0, u3_TH1=0):
    S, SR, C, CR, M1, M2, TH1, TH2, TC, Treg, IL10, IFNgamma, IL2 = y
    
    betaM2, betaTc, betaTh1CK2, betaTh1CK3, betaTh2, betaTreg, gammaC, gammaCR, gammaM1, gammaM2, gammaS, gammaTc, gammaTh1, gammaTh2, gammaTreg, deltaC, deltaCk1, deltaCk2, deltaCk3, deltaCR, deltaM1, deltaM2, deltaS, deltaTc, deltaTh1, deltaTh2, deltaTreg, lambdaM1, lambdaM2, lambdaTc1, lambdaTc2, lambdaTc3, lambdaTc4, lambdaTh1, lambdaTh2, lambdaTreg2, muC1, muC2, muS, muSR, muTcS, muTcTreg, muTh1Ck1, muTh1Ck3, muTregCk1, Cmax, CRmax, k1, k11, k2, k3, k4, k5, k6, k8, k9, ktc1, ktc2, ktc3, ktc4, mC, mS, p1, p2, r1, r2, tck, muM1Ck2, muM2Ck1, k7, k10 = params
    
    dS = (gammaS*((1-mS)*(1-p1-p2)))*S - (deltaS+(p2*gammaS)+gammaS*(mS*p1/2))*S - ((muS*S*IFNgamma)/(IFNgamma+k1)) - ((tck*S*TC)/(ktc1+TC)) - u2_S*S
    dSR = (gammaS*(1-p1-p2) - (deltaS+(p2*gammaS)))*SR + mS*gammaS*(1-p1/2-p2)*S - ((muSR*SR*IFNgamma)/(k2+IFNgamma)) - ((tck*SR*TC)/(ktc2+TC))
    dC = gammaC*(1-mC)*np.log((Cmax+log_eps)/(C+r1+log_eps))*C + gammaS*(p1+p2)*S - deltaC*C - mC*gammaC*C + (muC1*C*IL10)/(IL10+k3) - (muC2*C*IFNgamma)/(IFNgamma+k4) - (tck*C*TC)/(ktc3+TC) - u1*C - u2_C*C
    dCR = gammaCR*CR*np.log((CRmax+log_eps)/(CR+r2+log_eps)) + gammaS*SR*(p1+p2) + mC*gammaCR*C - deltaCR*CR + (muC1*CR*IL10)/(IL10+k5) - (muC2*CR*IFNgamma)/(IFNgamma+k6) - (tck*CR*TC)/(ktc4+TC) - u1*CR
    dM1 = gammaM1*M1*((C+CR)/(M1+lambdaM1)) - deltaM1*M1 + ((muM1Ck2*M1*IFNgamma)/(IFNgamma+k7))
    dM2 = gammaM2*M2*((C+CR)/(M2+lambdaM2)) - deltaM2*M2 + ((muM2Ck1*M2*IL10)/(IL10+k10))
    dTH1 = gammaTh1*((TH1*M1)/(lambdaTh1+TH1)) - deltaTh1*TH1 - ((muTh1Ck1*IL10*TH1)/(IL10+k8)) + ((muTh1Ck3*IL2*TH1)/(IL2+k9)) + u3_TH1
    dTH2 = gammaTh2*((TH2*M2)/(lambdaTh2+TH2)) - deltaTh2*TH2
    dTC = gammaTc*TC*((C+CR)/(TC+lambdaTc1)) + gammaTc*((TC*TH1)/(TC+lambdaTc4)) - muTcS*TC*((S+SR)/(TC+lambdaTc2)) - deltaTc*TC - muTcTreg*TC*((Treg)/(lambdaTc3+Treg)) + u3_Tc
    dTreg = gammaTreg*((Treg*M2)/(Treg+lambdaTreg2)) - deltaTreg*Treg + ((muTregCk1*IL10*Treg)/(Treg+k11))
    dIL10 = betaM2*M2 - deltaCk1*IL10 + betaTreg*Treg + betaTh2*TH2
    dIFNgamma = betaTh1CK2*TH1 + betaTc*TC - deltaCk2*IFNgamma
    dIL2 = betaTh1CK3*TH1 - deltaCk3*IL2
    
    return np.array([dS, dSR, dC, dCR, dM1, dM2, dTH1, dTH2, dTC, dTreg, dIL10, dIFNgamma, dIL2])

def event_equil(t, y, params):
    return np.linalg.norm(y) * 0.01

event_equil.terminal = True
event_equil.direction = -1

from scipy.integrate import odeint

def run_simulation(initial_conditions, params, t_final=1000, n_points=1000, u1=0, u2_S=0, u2_C=0, u3_Tc=0, u3_TH1=0):
    t_eval = np.linspace(0, t_final, n_points)
    
    def rhs_wrapper(y, t):
        return rhs(t, y, params, u1, u2_S, u2_C, u3_Tc, u3_TH1)

    y_eval = odeint(rhs_wrapper, initial_conditions, t_eval, 
                rtol=1e-4, atol=1e-6, mxstep=10000)  # Allow more steps
    
    state_names = ['S', 'SR', 'C', 'CR', 'M1', 'M2', 'TH1', 'TH2', 'TC', 'Treg', 'IL10', 'IFNgamma', 'IL2']
    df = pd.DataFrame(y_eval, columns=state_names)
    df.insert(0, 'time', t_eval)
    
    return df