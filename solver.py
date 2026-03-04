import sys
import shutil
import os
import numpy as np
from scipy.integrate import solve_ivp
import pandas as pd
from model import rhs, get_default_ic, get_params

# Inputs (customize here or via CLI)
t_final = 'equil' # 800.0  
# Or 'equil' for auto-detect

output_file = 'simulation_results.csv'
# ranges = None  # e.g., {0: (0.5, 1.5), 1: (0.8, 1.2)} for param indices
initial_condns = get_default_ic()  # Or custom np.array

params = get_params()

def event_equil(t, y, params): 
    return np.linalg.norm(y) * 0.01  # Stop near zero activity

event_equil.terminal = True
event_equil.direction = -1

sol = solve_ivp(rhs, 
                [0, t_final], 
                initial_condns, 
                args=(params,), 
                method='BDF', rtol=1e-6, atol=1e-8, 
                events=event_equil if t_final == 'equil' else None, 
                dense_output=True)

if not sol.success:
    print("Solver failed:", sol.message)
    sys.exit(1)

t_eval = np.linspace(0, sol.t[-1], 1000)
y_eval = sol.sol(t_eval).T

# States columns
state_names = ['S', 'SR', 'C', 'CR', 'M1', 'M2', 'TH1', 'TH2', 'TC', 'Treg', 'IL10', 'IFNgamma', 'IL2']
df = pd.DataFrame(y_eval, columns=state_names)
df.insert(0, 'time', t_eval)

df.to_csv(output_file, index=False)
print(f"Results saved to {output_file} (t_end={sol.t[-1]:.2f})")

# deleting the __pycache folder
if os.path.exists('__pycache__'):
    shutil.rmtree('__pycache__')