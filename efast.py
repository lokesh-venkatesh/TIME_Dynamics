import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.integrate import odeint
from scipy import stats
from multiprocessing import Pool
import functools
from model import get_params, get_default_ic, rhs, PARAM_NAMES

# Configuration matching paper
N_PARAMS = 71
N_SAMPLES = 10  # samples per search curve
N_RESAMPLES = 3  # 3 times resampling
TOTAL_SIMS = (N_PARAMS + 1) * N_SAMPLES * N_RESAMPLES  # 21,600 simulations
TIME_POINTS = [2, 4, 6, 10, 100, 600]  # days - tumor development stages
PARAM_LIST = sorted(PARAM_NAMES.keys(), key=lambda x: PARAM_NAMES[x])
SPECIES_INDICES = {'S': 0, 'SR': 1, 'C': 2, 'CR': 3}
N_CORES = 8  # Use 8 out of 12 cores

# ============================================================================
# Standalone worker function for multiprocessing
# ============================================================================
def analyze_parameter_worker(param_idx, omega_set, t_point, species_idx, n_samples):
    """
    Standalone worker function for parameter analysis.
    Must be module-level and pickleable for multiprocessing.
    """
    sensitivities = []
    
    for resample_idx, omega in enumerate(omega_set):
        try:
            # Generate base curve
            s = np.linspace(0, 2 * np.pi, n_samples)
            samples_base = np.zeros((n_samples, N_PARAMS))
            
            for i in range(N_PARAMS):
                curve = np.sin(omega[i] * s)
                samples_base[:, i] = (curve + 1.0) / 2.0
            
            # Transform to parameters
            params_base = transform_to_parameters(samples_base)
            
            # Run base model
            Y_base = []
            for j in range(n_samples):
                y_val = run_model(params_base[j], t_point, species_idx)
                Y_base.append(y_val)
            Y_base = np.array(Y_base)
            
            # Generate perturbed curve
            samples_perturb = samples_base.copy()
            samples_perturb[:, param_idx] = (samples_perturb[:, param_idx] + 0.2) % 1.0
            params_perturb = transform_to_parameters(samples_perturb)
            
            # Run perturbed model
            Y_perturb = []
            for j in range(n_samples):
                y_val = run_model(params_perturb[j], t_point, species_idx)
                Y_perturb.append(y_val)
            Y_perturb = np.array(Y_perturb)
            
            # Compute sensitivity
            S1, ST = compute_variance_indices(Y_base, Y_perturb)
            if not np.isnan(S1):
                sensitivities.append(S1)
        
        except Exception as e:
            continue
    
    if len(sensitivities) == 0:
        return np.nan, 1.0
    
    mean_S1 = np.mean(sensitivities)
    pval = compute_pvalue(sensitivities)
    
    return mean_S1, pval

# ============================================================================
# Helper functions (module-level for pickling)
# ============================================================================
def run_model(params, t_point, species_idx):
    """Run ODE model and return species value at time t"""
    t_eval = np.linspace(0, t_point, max(50, int(t_point * 2)))
    
    def rhs_wrapper(y, t):
        return rhs(t, y, params)
    
    ic = get_default_ic()
    try:
        y_eval = odeint(rhs_wrapper, ic, t_eval, 
                       rtol=1e-4, atol=1e-6, mxstep=10000, full_output=False)
        return y_eval[-1, species_idx]
    except:
        return np.nan

def transform_to_parameters(normalized_samples, param_ranges=None):
    """Transform normalized samples [0,1] to actual parameter values"""
    base_params = get_params()
    
    if param_ranges is None:
        param_ranges = {}
        for i in range(N_PARAMS):
            val = base_params[i]
            if val == 0:
                param_ranges[i] = (1e-16, 1e-14)
            else:
                param_ranges[i] = (val * 0.5, val * 1.5)
    
    param_samples = np.zeros_like(normalized_samples)
    for i in range(N_PARAMS):
        low, high = param_ranges[i]
        param_samples[:, i] = low + normalized_samples[:, i] * (high - low)
    
    return param_samples

def compute_variance_indices(Y_base, Y_perturb):
    """Compute first-order and total sensitivity indices"""
    if np.any(np.isnan(Y_base)) or np.any(np.isnan(Y_perturb)):
        return np.nan, np.nan
    
    f_0 = np.mean(Y_base)
    D_base = np.var(Y_base)
    D_perturb = np.var(Y_perturb)
    
    if D_base < 1e-15:
        return 0.0, 0.0
    
    S1 = (D_perturb - D_base) / D_base if D_base > 0 else 0.0
    S1 = np.clip(abs(S1), 0, 1)
    
    ST = 1.0 - (D_perturb / D_base) if D_base > 0 else 0.0
    ST = np.clip(ST, 0, 1)
    
    return S1, ST

def compute_pvalue(sensitivities):
    """Compute p-value from sensitivity distribution"""
    if len(sensitivities) < 2 or np.all(np.isnan(sensitivities)):
        return 1.0
    
    sensitivities = sensitivities[~np.isnan(sensitivities)]
    if len(sensitivities) < 2:
        return 1.0
    
    mean_s = np.mean(sensitivities)
    std_s = np.std(sensitivities)
    
    if std_s > 1e-10:
        t_stat = mean_s / (std_s / np.sqrt(len(sensitivities)))
        pval = 1 - stats.norm.cdf(abs(t_stat))
    else:
        pval = 1.0 if mean_s < 0.01 else 0.0
    
    return pval

# ============================================================================
# Original eFASTAnalyzer class (unchanged)
# ============================================================================
class eFASTAnalyzer:
    """Extended Fourier Amplitude Sensitivity Test"""
    
    def __init__(self, n_samples=100, n_resamples=5):
        self.n_samples = n_samples
        self.n_resamples = n_resamples
        self.results = {}
    
    def generate_omega_set(self):
        """Generate random frequency set for each resample"""
        omega_set = []
        for _ in range(self.n_resamples):
            omega = np.random.choice(range(1, N_PARAMS + 10), size=N_PARAMS, replace=False)[:N_PARAMS]
            omega_set.append(omega)
        return omega_set
    
    def generate_samples(self, omega):
        """Generate eFAST curve samples using sine transformation"""
        s = np.linspace(0, 2 * np.pi, self.n_samples)
        samples = np.zeros((self.n_samples, N_PARAMS))
        
        for i in range(N_PARAMS):
            curve = np.sin(omega[i] * s)
            samples[:, i] = (curve + 1.0) / 2.0
        
        return samples
    
    def transform_to_parameters(self, normalized_samples, param_ranges=None):
        """Transform normalized samples [0,1] to actual parameter values"""
        base_params = get_params()
        
        if param_ranges is None:
            param_ranges = {}
            for i in range(N_PARAMS):
                val = base_params[i]
                if val == 0:
                    param_ranges[i] = (1e-16, 1e-14)
                else:
                    param_ranges[i] = (val * 0.5, val * 1.5)
        
        param_samples = np.zeros_like(normalized_samples)
        for i in range(N_PARAMS):
            low, high = param_ranges[i]
            param_samples[:, i] = low + normalized_samples[:, i] * (high - low)
        
        return param_samples
    
    def run_model(self, params, t_point, species_idx):
        """Run ODE model and return species value at time t"""
        t_eval = np.linspace(0, t_point, max(50, int(t_point * 2)))
        
        def rhs_wrapper(y, t):
            return rhs(t, y, params)
        
        ic = get_default_ic()
        try:
            y_eval = odeint(rhs_wrapper, ic, t_eval, 
                           rtol=1e-4, atol=1e-6, mxstep=10000, full_output=False)
            return y_eval[-1, species_idx]
        except:
            return np.nan
    
    def compute_variance_indices(self, Y_base, Y_perturb):
        """Compute first-order and total sensitivity indices"""
        if np.any(np.isnan(Y_base)) or np.any(np.isnan(Y_perturb)):
            return np.nan, np.nan
        
        f_0 = np.mean(Y_base)
        D_base = np.var(Y_base)
        D_perturb = np.var(Y_perturb)
        
        if D_base < 1e-15:
            return 0.0, 0.0
        
        S1 = (D_perturb - D_base) / D_base if D_base > 0 else 0.0
        S1 = np.clip(abs(S1), 0, 1)
        
        ST = 1.0 - (D_perturb / D_base) if D_base > 0 else 0.0
        ST = np.clip(ST, 0, 1)
        
        return S1, ST
    
    def compute_pvalue(self, sensitivities):
        """Compute p-value from sensitivity distribution"""
        if len(sensitivities) < 2 or np.all(np.isnan(sensitivities)):
            return 1.0
        
        sensitivities = sensitivities[~np.isnan(sensitivities)]
        if len(sensitivities) < 2:
            return 1.0
        
        mean_s = np.mean(sensitivities)
        std_s = np.std(sensitivities)
        
        if std_s > 1e-10:
            t_stat = mean_s / (std_s / np.sqrt(len(sensitivities)))
            pval = 1 - stats.norm.cdf(abs(t_stat))
        else:
            pval = 1.0 if mean_s < 0.01 else 0.0
        
        return pval
    
    def filter_results(self, p_threshold=0.05):
        """Filter significant parameters (p < threshold)"""
        filtered = {}
        
        for species_name in SPECIES_INDICES.keys():
            filtered[species_name] = {}
            
            for t_point in TIME_POINTS:
                data = self.results[species_name][t_point]
                df = pd.DataFrame({
                    'Parameter': data['params'],
                    'S1': data['S1'],
                    'pval': data['pval']
                })
                
                df_sig = df[df['pval'] < p_threshold].sort_values('S1', ascending=False)
                filtered[species_name][t_point] = df_sig
        
        return filtered
    
    def plot_heatmaps(self, filtered, save_dir='data/sensitivity'):
        """Create 2D seaborn heatmaps"""
        os.makedirs(save_dir, exist_ok=True)
        
        species_info = {
            'S': 'Stem Cells',
            'SR': 'Stem Resistant Cells',
            'C': 'Cancer Cells',
            'CR': 'Cancer Resistant Cells'
        }
        
        for species_name, species_title in species_info.items():
            print(f"\nCreating heatmap for {species_title}...")
            
            all_params = sorted(set().union(*[
                set(filtered[species_name][tp]['Parameter'].tolist()) 
                for tp in TIME_POINTS
            ]))
            
            if len(all_params) == 0:
                print(f"  No significant parameters found for {species_name}")
                continue
            
            heatmap_data = np.zeros((len(TIME_POINTS), len(all_params)))
            heatmap_data[:] = np.nan
            
            for i, t_point in enumerate(TIME_POINTS):
                df_t = filtered[species_name][t_point]
                for j, param in enumerate(all_params):
                    s1_val = df_t[df_t['Parameter'] == param]['S1'].values
                    if len(s1_val) > 0:
                        heatmap_data[i, j] = s1_val[0]
            
            heatmap_df = pd.DataFrame(
                heatmap_data,
                index=[f'{int(tp)} days' for tp in TIME_POINTS],
                columns=all_params
            )
            
            fig, ax = plt.subplots(figsize=(18, 6))
            
            sns.heatmap(heatmap_df, 
                       cmap='YlOrRd', 
                       cbar_kws={'label': 'Sensitivity Index (S1)'},
                       ax=ax,
                       linewidths=0.5,
                       linecolor='gray')
            
            ax.set_title(f'eFAST Parameter Sensitivity: {species_title}\n(p < 0.05)', 
                        fontsize=14, fontweight='bold', pad=20)
            ax.set_xlabel('Parameter', fontsize=12, fontweight='bold')
            ax.set_ylabel('Time Point', fontsize=12, fontweight='bold')
            
            plt.xticks(rotation=45, ha='right', fontsize=10)
            plt.yticks(fontsize=11)
            plt.tight_layout()
            
            plt.savefig(f'{save_dir}/sensitivity_heatmap_{species_name}.png', 
                       dpi=300, bbox_inches='tight')
            print(f"  ✓ Saved: sensitivity_heatmap_{species_name}.png")
            plt.close()
    
    def save_results_csv(self, filtered, save_dir='data/sensitivity'):
        """Save results to CSV files"""
        os.makedirs(save_dir, exist_ok=True)
        
        for species_name in SPECIES_INDICES.keys():
            all_results = []
            
            for t_point in TIME_POINTS:
                df = filtered[species_name][t_point].copy()
                df['Time_days'] = t_point
                all_results.append(df)
            
            if all_results:
                combined_df = pd.concat(all_results, ignore_index=True)
                combined_df = combined_df[['Time_days', 'Parameter', 'S1', 'pval']]
                combined_df = combined_df.sort_values(['Time_days', 'S1'], ascending=[True, False])
                
                filename = f'{save_dir}/efast_results_{species_name}.csv'
                combined_df.to_csv(filename, index=False)
                print(f"✓ Saved: {filename}")

# ============================================================================
# Parallelized Analysis Method
# ============================================================================
class ParalleleFASTAnalyzer(eFASTAnalyzer):
    """Extended Fourier Amplitude Sensitivity Test with Parallelization"""
    
    def run_analysis(self):
        """Run full sensitivity analysis with parallelization across parameters"""
        print("\n" + "="*80)
        print("eFAST PARAMETER SENSITIVITY ANALYSIS (PARALLELIZED)")
        print(f"Total simulations: {TOTAL_SIMS:,}")
        print(f"Parameters: {N_PARAMS}, Samples/curve: {N_SAMPLES}, Resamples: {N_RESAMPLES}")
        print(f"Cores used: {N_CORES}")
        print(f"Timepoints: {TIME_POINTS}")
        print("="*80)
        
        omega_set = self.generate_omega_set()
        print(f"\n✓ Generated {self.n_resamples} frequency sets\n")
        
        for species_name, species_idx in SPECIES_INDICES.items():
            print(f"{'='*80}")
            print(f"ANALYZING: {species_name} (Species Index {species_idx})")
            print(f"{'='*80}")
            
            self.results[species_name] = {}
            
            for t_idx, t_point in enumerate(TIME_POINTS):
                print(f"\n  Time = {t_point} days ({t_idx+1}/{len(TIME_POINTS)})")
                print(f"  Parallelizing across {N_PARAMS} parameters on {N_CORES} cores...")
                
                # Create partial function with fixed arguments
                worker_fn = functools.partial(
                    analyze_parameter_worker,
                    omega_set=omega_set,
                    t_point=t_point,
                    species_idx=species_idx,
                    n_samples=self.n_samples
                )
                
                # Run parameter analysis in parallel
                with Pool(processes=N_CORES) as pool:
                    results = pool.map(worker_fn, range(N_PARAMS))
                
                # Unpack results
                s1_values, pval_values = zip(*results)
                
                # Store results
                self.results[species_name][t_point] = {
                    'S1': np.array(s1_values),
                    'pval': np.array(pval_values),
                    'params': PARAM_LIST
                }
                
                # Summary
                pval_array = np.array(pval_values)
                n_sig = np.sum(pval_array < 0.05)
                sig_S1 = [s for s, p in zip(s1_values, pval_values) if p < 0.05]
                mean_S1 = np.mean(sig_S1) if sig_S1 else 0.0
                
                print(f"  ✓ Found {n_sig} significant parameters (p < 0.05)")
                if sig_S1:
                    print(f"    Mean S1 for significant params: {mean_S1:.4f}")

def main():
    analyzer = ParalleleFASTAnalyzer(n_samples=N_SAMPLES, n_resamples=N_RESAMPLES)
    
    # Run full sensitivity analysis with parallelization
    analyzer.run_analysis()
    
    # Filter significant parameters
    filtered = analyzer.filter_results(p_threshold=0.05)
    
    # Create heatmap visualizations
    analyzer.plot_heatmaps(filtered)
    
    # Save results
    print("\n" + "="*80)
    print("SAVING RESULTS")
    print("="*80)
    analyzer.save_results_csv(filtered)
    
    print("\n" + "="*80)
    print("SENSITIVITY ANALYSIS COMPLETE")
    print("="*80 + "\n")

if __name__ == '__main__':
    # Required for Windows multiprocessing
    from multiprocessing import freeze_support
    freeze_support()
    
    main()