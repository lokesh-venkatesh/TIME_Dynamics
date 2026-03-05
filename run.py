import os
import shutil
import numpy as np
import pandas as pd
from model import get_params, get_default_ic, run_simulation

os.makedirs('data', exist_ok=True)

def run_default():
    params = get_params()
    ic = get_default_ic()
    return run_simulation(ic, params, t_final=10000) # 'equil'

def run_modified(param_indices, param_values, t_final=10000):
    params = get_params()
    for idx, val in zip(param_indices, param_values):
        print(params[idx])
        params[idx] = val
    ic = get_default_ic()
    return run_simulation(ic, params, t_final=t_final)

def run_sweep(param_sweep_dict, output_prefix=None):
    param_indices = list(param_sweep_dict.keys())
    param_ranges = [param_sweep_dict[idx] for idx in param_indices]
    
    if len(param_indices) == 1:
        idx = param_indices[0]
        values = param_ranges[0]
        results = {}
        for val in values:
            try:
                params = get_params()
                params[idx] = val
                ic = get_default_ic()
                results[val] = run_simulation(ic, params, t_final=10000)
                print(f"✓ param[{idx}] = {val:.6f}")
            except RuntimeError as e:
                print(f"✗ param[{idx}] = {val:.6f} - {e}")
        
        if output_prefix:
            for val, df in results.items():
                if df is not None:
                    df.to_csv(f"data/{output_prefix}_param{idx}_{val:.6f}.csv", index=False)
        return results
    
    else:
        results = {}
        for combo in np.ndindex(tuple(len(r) for r in param_ranges)):
            params = get_params()
            param_combo = {}
            for i, idx in enumerate(param_indices):
                val = param_ranges[i][combo[i]]
                params[idx] = val
                param_combo[idx] = val
            
            ic = get_default_ic()
            key = tuple((param_indices[i], param_ranges[i][combo[i]]) for i in range(len(param_indices)))
            results[key] = run_simulation(ic, params, t_final=10000)
            
            if output_prefix:
                filename = output_prefix
                for idx, val in key:
                    filename += f"_p{idx}_{val:.6e}"
                results[key].to_csv(f"data/{filename}.csv", index=False)
        
        return results

if __name__ == '__main__':


    # Example usage:
    
    # 1. Run with default parameters
    # df_default = run_default()
    # df_default.to_csv('simulation_results.csv', index=False)
    # print(f"Default run completed, saved to results_default.csv")
    
    # 2. Modify specific parameters
    # df_mod = run_modified([0, 1], [2e-15, 2e-8])
    # df_mod = run_modified([2], [0.1], t_final=1000)
    # df_mod.to_csv('results_for_figure_2i.csv', index=False)
    
    # 3. Parameter sweep (single parameter)
    # results = run_sweep({62: np.linspace(0.2, 0.4, 100)}, output_prefix='sweep_param62')
    # print(results)
    
    # 4. Parameter sweep (multiple parameters)
    # results = run_sweep({0: np.linspace(1e-16, 1e-14, 3), 
    #                      1: np.linspace(1e-9, 1e-7, 3)}, 
    #                     output_prefix='sweep_multi')

    # deleting the __pycache folder
    if os.path.exists('__pycache__'):
        shutil.rmtree('__pycache__')