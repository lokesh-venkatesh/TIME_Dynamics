import sys
import numpy as np
from scipy.integrate import solve_ivp
import pandas as pd
from model import rhs
from params import get_params, get_default_ic

# Inputs (customize here or via CLI)
t_final = 800.0  # Or 'equil' for auto-detect
output_file = 'simulation_results.csv'
ranges = None  # e.g., {0: (0.5, 1.5), 1: (0.8, 1.2)} for param indices
ic = get_default_ic()  # Or custom np.array

params = get_params(ranges)

def event_equil(t, y, params): return np.linalg.norm(y) * 0.01  # Stop near zero activity
event_equil.terminal = True
event_equil.direction = -1

sol = solve_ivp(rhs, [0, t_final], ic, args=(params,), method='BDF', rtol=1e-6, atol=1e-8, events=event_equil if t_final == 'equil' else None, dense_output=True)

if not sol.success:
    print("Solver failed:", sol.message)
    sys.exit(1)

t_eval = np.linspace(0, sol.t[-1], 1000)
y_eval = sol.sol(t_eval).T

# States columns
state_names = ['S', 'SR', 'C', 'CR', 'M1', 'M2', 'TH1', 'TH2', 'TC', 'Treg', 'IL10', 'IFNgamma', 'IL2']
df = pd.DataFrame(y_eval, columns=state_names)
df.insert(0, 'time', t_eval)

# Append params and ICs
param_names = ['betaM2', 'betaTc', 'betaTh1CK2', 'betaTh1CK3', 'betaTh2', 'betaTreg', 'gammaC', 'gammaCR', 'gammaM1', 'gammaM2', 'gammaS', 'gammaTc', 'gammaTh1', 'gammaTh2', 'gammaTreg',
               'deltaC', 'deltaCk1', 'deltaCk2', 'deltaCk3', 'deltaCR', 'deltaM1', 'deltaM2', 'deltaS', 'deltaTc', 'deltaTh1', 'deltaTh2', 'deltaTreg',
               'lambdaM1', 'lambdaM2', 'lambdaTc1', 'lambdaTc2', 'lambdaTc3', 'lambdaTc4', 'lambdaTh1', 'lambdaTh2', 'lambdaTreg2',
               'muC1', 'muC2', 'muS', 'muSR', 'muTcS', 'muTcTreg', 'muTh1Ck1', 'muTh1Ck3', 'muTregCk1',
               'Cmax', 'CRmax',
               'k1', 'k11', 'k2', 'k3', 'k4', 'k5', 'k6', 'k8', 'k9', 'ktc1', 'ktc2', 'ktc3', 'ktc4', 'mC', 'mS', 'p1', 'p2', 'r1', 'r2', 'tck', 'muM1Ck2', 'muM2Ck1', 'k7', 'k10']

# COMMENTING OUT FOR THE TIME-BEING, SINCE I DON'T WANT PARAM INFORMATION ABHI 
#for i, name in enumerate(param_names):
#    df[name] = params[i]

ic_names = [f'ic_{name}' for name in state_names]
for i, name in enumerate(ic_names):
    df[name] = ic[i]

df.to_csv(output_file, index=False)
print(f"Results saved to {output_file} (t_end={sol.t[-1]:.2f})")
