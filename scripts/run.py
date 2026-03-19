import os
import shutil
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from model import get_params, get_default_ic, run_simulation, get_param_index 

def run_default(T='equil'):
    params = get_params()
    ic = get_default_ic()
    return run_simulation(ic, params, t_final=T) 

def run_modified(param_keys, param_values, t_final=800): # Accept parameter names or indices
    params = get_params()
    for key, val in zip(param_keys, param_values): # Convert names to indices if needed
        idx = get_param_index(key)
        params[idx] = val
    ic = get_default_ic()
    return run_simulation(ic, params, t_final=t_final)

def run_sweep(param_sweep_dict, output_prefix=None): # Accept parameter names or indices as keys
    param_keys = list(param_sweep_dict.keys())
    param_indices = [get_param_index(key) for key in param_keys] # Convert all keys to indices
    param_ranges = [param_sweep_dict[key] for key in param_keys]
    
    if len(param_keys) == 1:
        key = param_keys[0]
        idx = param_indices[0]
        values = param_ranges[0]
        results = {}
        for val in values:
            try:
                params = get_params()
                params[idx] = val
                ic = get_default_ic()
                results[val] = run_simulation(ic, params, t_final=2000)
                print(f"✓ {key if isinstance(key, str) else 'param['+str(key)+']'} = {val:.6f}") # Print name if available
            except RuntimeError as e:
                print(f"✗ {key if isinstance(key, str) else 'param['+str(key)+']'} = {val:.6f} - {e}")
        
        if output_prefix:
            for val, df in results.items():
                if df is not None:
                    ...
                    # df.to_csv(f"data/more_data/{output_prefix}_param{idx}_{val:.6f}.csv", index=False)
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
            print(f"Running combo: " + ", ".join([f"{key if isinstance(key, str) else 'param['+str(key)+']'}={val:.6f}" for key, val in zip(param_keys, [param_combo[idx] for idx in param_indices])]))
            
            ic = get_default_ic()
            key = tuple((param_indices[i], param_ranges[i][combo[i]]) for i in range(len(param_indices)))
            results[key] = run_simulation(ic, params, t_final=2000)
            
            if output_prefix:
                filename = output_prefix
                for idx, val in key:
                    filename += f"_p{idx}_{val:.6e}"
                # results[key].to_csv(f"data/more_data/{filename}.csv", index=False)
        
        return results

if __name__ == '__main__':
    os.makedirs('data', exist_ok=True)

    # Example usage:
    
    # 1. Run with default parameters
    df = run_default()
    plt.figure(figsize=(7,5))


    plt.plot(df["time"], df["S"], label="S", linewidth=2)
    plt.plot(df["time"], df["SR"], label="SR", linewidth=2)
    plt.plot(df["time"], df["C"], label="C", linewidth=2)
    plt.plot(df["time"], df["CR"], label="CR", linewidth=2)
    plt.plot(df["time"], df["M1"], label="M1", linewidth=2)
    plt.plot(df["time"], df["M2"], label="M2", linewidth=2)
    plt.plot(df["time"], df["TH1"], label="TH1", linewidth=2)
    plt.plot(df["time"], df["TH2"], label="TH2", linewidth=2)
    plt.plot(df["time"], df["TC"], label="TC", linewidth=2)
    plt.plot(df["time"], df["Treg"], label="Treg", linewidth=2)
    
    plt.plot(df["time"], df["IL10"], color="red", linewidth=2, label="IL10")
    plt.plot(df["time"], df["IFNgamma"], color="green", linewidth=2, label="IFN-γ")
    plt.plot(df["time"], df["IL2"], color="black", linewidth=2, label="IL2")

    plt.xlabel("Time (days)", fontsize=12)
    plt.ylabel("Cell Density (cells/ml)", fontsize=12)
    plt.legend(ncol=2, fontsize=8)
    plt.xlim(left=0)
    plt.ylim(bottom=0)

    plt.tight_layout()
    # plt.savefig("figure1/figure_1c.png", dpi=300)
    # plt.close()
    plt.show()

    # df_default.to_csv('simulation_results.csv', index=False)
    # print(f"Default run completed, saved to results_default.csv")
    
    # 2. Modify specific parameters BY NAME
    # df_mod = run_modified(['gammaM1', 'gammaM2'], [0.3, 0.1])
    # df_mod = run_modified(['gammaM1'], [0.25], t_final=1000)
    # df_mod.to_csv('results_for_figure_2i.csv', index=False)
    
    # 3. Parameter sweep by name (single parameter)
    # results = run_sweep({'gammaM1': np.linspace(0.2, 0.4, 100)}, output_prefix='sweep_gammaM1')
    # print(results)
    
    # 4. Parameter sweep (multiple parameters)
    # results = run_sweep({'betaM2': np.linspace(1e-16, 1e-14, 3), 
    #                      'betaTc': np.linspace(1e-9, 1e-7, 3)}, 
    #                     output_prefix='sweep_multi')

    if os.path.exists('__pycache__'):
        shutil.rmtree('__pycache__')