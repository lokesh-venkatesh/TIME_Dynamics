import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.stats import norm
from scipy.spatial.distance import euclidean
import seaborn as sns
from model import get_params, get_default_ic, rhs, PARAM_NAMES, get_param_index

# ============================================================================
# Experimental Data: SGC7901 Gastric Cancer Cell Line (7 days)
# ============================================================================
# Replace with your actual experimental data if available
# Format: time (days), cancer cell count (cells)
EXPERIMENTAL_DATA = {
    'time': np.array([0, 1, 2, 3, 4, 5, 6, 7]),
    'cancer_cells': np.array([1.0, 2.1, 4.8, 10.5, 23.0, 50.0, 108.0, 235.0])  # Relative counts
}

# Unknown parameters to estimate (indices in parameter array)
# These are parameters not curated from literature and need estimation
UNKNOWN_PARAMS = [
    ('gammaC', 6),      # Cancer cell growth rate
    ('deltaC', 15),     # Cancer cell death rate
    ('muC1', 36),       # IL10-mediated cancer suppression
    ('muC2', 37),       # IFNgamma-mediated cancer suppression
    ('Cmax', 45),       # Cancer cell carrying capacity
]

# Prior distributions (mean, std for normal prior)
PRIORS = {
    'gammaC': (0.1, 0.05),
    'deltaC': (0.8, 0.4),
    'muC1': (0.75, 0.3),
    'muC2': (0.9, 0.4),
    'Cmax': (1e10, 5e9),
}

# Biological feasible ranges
PARAM_BOUNDS = {
    'gammaC': (0.01, 0.5),
    'deltaC': (0.01, 2.0),
    'muC1': (0.1, 1.5),
    'muC2': (0.1, 1.5),
    'Cmax': (1e9, 1e11),
}

MCMC_CONFIG = {
    'n_iterations': 50000,      # 5 lakh = 500,000 in paper; reduced for demo
    'n_burn': 10000,            # Burn-in period
    'n_chains': 4,              # Multiple chains for convergence diagnostics
    'dr_levels': 2,             # Delayed rejection levels
    'adaptation_interval': 100, # Adapt covariance every N iterations
}

# ============================================================================
# Likelihood Function
# ============================================================================
class CancerCellDataLikelihood:
    """Likelihood function for fitting cancer cell growth to experimental data"""
    
    def __init__(self, exp_data, unknown_params_dict):
        self.exp_time = exp_data['time']
        self.exp_cells = exp_data['cancer_cells']
        self.unknown_params_dict = unknown_params_dict
        self.measurement_error = 0.1  # 10% measurement error
    
    def simulate_cancer_cells(self, params):
        """Simulate cancer cell dynamics"""
        ic = get_default_ic()
        
        def rhs_wrapper(y, t):
            return rhs(t, y, params)
        
        try:
            y_eval = odeint(rhs_wrapper, ic, self.exp_time,
                          rtol=1e-4, atol=1e-6, mxstep=10000)
            # Return C (cancer cells) column
            return y_eval[:, 2]  # Column 2 is 'C'
        except:
            return np.full_like(self.exp_time, np.inf)
    
    def log_likelihood(self, sim_cells):
        """Compute log-likelihood: -0.5 * sum((y_sim - y_obs)^2 / error^2)"""
        if np.any(np.isnan(sim_cells)) or np.any(sim_cells < 0):
            return -np.inf
        
        # Normalize both to same scale for comparison
        exp_norm = self.exp_cells / (np.max(self.exp_cells) + 1e-10)
        sim_norm = sim_cells / (np.max(sim_cells) + 1e-10)
        
        error_var = (self.measurement_error ** 2) * np.mean(exp_norm ** 2)
        residuals = (sim_norm - exp_norm) ** 2
        
        log_like = -0.5 * np.sum(residuals / error_var)
        return log_like
    
    def evaluate(self, param_vector):
        """Evaluate likelihood for given parameter vector"""
        base_params = get_params()
        
        # Update unknown parameters
        for i, param_name in enumerate(list(self.unknown_params_dict.keys())):
            param_idx = self.unknown_params_dict[param_name]
            base_params[param_idx] = param_vector[i]
        
        # Simulate and compute likelihood
        sim_cells = self.simulate_cancer_cells(base_params)
        log_like = self.log_likelihood(sim_cells)
        
        return log_like

# ============================================================================
# MCMC-DRAM Sampler
# ============================================================================
class MCMCDRAMSampler:
    """Markov Chain Monte Carlo with Delayed Rejection Adaptive Metropolis"""
    
    def __init__(self, likelihood, unknown_params_dict, priors, bounds):
        self.likelihood = likelihood
        self.unknown_params_dict = unknown_params_dict
        self.priors = priors
        self.bounds = bounds
        self.param_names = list(unknown_params_dict.keys())
        self.n_params = len(self.param_names)
        
        # Initialize covariance matrix
        self.cov = np.eye(self.n_params) * 0.1
        self.cov_adapt_counter = 0
    
    def log_prior(self, param_vector):
        """Evaluate log prior (product of normal distributions)"""
        log_p = 0.0
        for i, param_name in enumerate(self.param_names):
            mean, std = self.priors[param_name]
            log_p += norm.logpdf(param_vector[i], loc=mean, scale=std)
        return log_p
    
    def is_within_bounds(self, param_vector):
        """Check if parameters are within biological feasible ranges"""
        for i, param_name in enumerate(self.param_names):
            low, high = self.bounds[param_name]
            if param_vector[i] < low or param_vector[i] > high:
                return False
        return True
    
    def log_posterior(self, param_vector):
        """Evaluate log posterior = log prior + log likelihood"""
        if not self.is_within_bounds(param_vector):
            return -np.inf
        
        log_p = self.log_prior(param_vector)
        log_like = self.likelihood.evaluate(param_vector)
        
        return log_p + log_like
    
    def metropolis_step(self, current_params, scale=1.0, dr_level=0):
        """Single Metropolis step with optional delayed rejection"""
        # Propose new parameters
        proposal = current_params + scale * np.random.multivariate_normal(
            np.zeros(self.n_params), self.cov
        )
        
        # Compute acceptance ratio
        log_alpha = self.log_posterior(proposal) - self.log_posterior(current_params)
        
        # Accept or reject
        if np.log(np.random.uniform()) < log_alpha:
            return proposal, True
        else:
            # Delayed rejection: try with smaller step
            if dr_level < MCMC_CONFIG['dr_levels']:
                smaller_scale = scale / 2.0
                return self.metropolis_step(current_params, smaller_scale, dr_level + 1)
            else:
                return current_params, False
    
    def run_chain(self, initial_params=None, chain_id=0):
        """Run single MCMC chain"""
        if initial_params is None:
            # Initialize from prior means
            initial_params = np.array([self.priors[pn][0] for pn in self.param_names])
        
        chain = np.zeros((MCMC_CONFIG['n_iterations'], self.n_params))
        log_post = np.zeros(MCMC_CONFIG['n_iterations'])
        acceptance_rate = 0
        
        current_params = initial_params.copy()
        current_log_post = self.log_posterior(current_params)
        
        print(f"  Chain {chain_id}: Running {MCMC_CONFIG['n_iterations']} iterations...")
        
        for it in range(MCMC_CONFIG['n_iterations']):
            # Metropolis step with DRAM
            current_params, accepted = self.metropolis_step(current_params)
            if accepted:
                acceptance_rate += 1
            
            current_log_post = self.log_posterior(current_params)
            
            chain[it] = current_params
            log_post[it] = current_log_post
            
            # Adaptive covariance update
            if (it + 1) % MCMC_CONFIG['adaptation_interval'] == 0:
                self.cov = np.cov(chain[:it+1].T)
                self.cov += 1e-6 * np.eye(self.n_params)  # Add small positive value
            
            if (it + 1) % 5000 == 0:
                acc_rate = acceptance_rate / (it + 1)
                print(f"    Iteration {it+1}: Acceptance rate = {acc_rate:.2%}, "
                      f"Log-posterior = {current_log_post:.2f}")
        
        final_acc_rate = acceptance_rate / MCMC_CONFIG['n_iterations']
        print(f"  ✓ Chain {chain_id} complete. Final acceptance rate: {final_acc_rate:.2%}\n")
        
        return chain, log_post
    
    def run_parallel_chains(self):
        """Run multiple chains in parallel for convergence diagnostics"""
        print("\n" + "="*80)
        print("RUNNING MCMC-DRAM PARAMETER ESTIMATION")
        print(f"Parameters: {self.param_names}")
        print(f"Iterations: {MCMC_CONFIG['n_iterations']:,}, Chains: {MCMC_CONFIG['n_chains']}")
        print("="*80 + "\n")
        
        chains_list = []
        log_post_list = []
        
        for chain_id in range(MCMC_CONFIG['n_chains']):
            chain, log_post = self.run_chain(chain_id=chain_id)
            chains_list.append(chain)
            log_post_list.append(log_post)
        
        return np.array(chains_list), np.array(log_post_list)
    
    def gelman_rubin_diagnostic(self, chains):
        """Compute Gelman-Rubin statistic for convergence assessment"""
        n_chains = chains.shape[0]
        n_samples = chains.shape[1]
        n_params = chains.shape[2]
        
        gr_stats = np.zeros(n_params)
        
        for p in range(n_params):
            # Between-chain variance
            chain_means = np.mean(chains[:, :, p], axis=1)
            B = n_samples * np.var(chain_means)
            
            # Within-chain variance
            chain_vars = np.var(chains[:, :, p], axis=1)
            W = np.mean(chain_vars)
            
            # Compute Gelman-Rubin statistic
            var_hat = ((n_samples - 1) / n_samples) * W + (1 / n_samples) * B
            gr_stats[p] = np.sqrt(var_hat / W)
        
        return gr_stats
    
    def plot_diagnostics(self, chains, log_post, save_dir='data/mcmc'):
        """Create convergence diagnostic plots"""
        os.makedirs(save_dir, exist_ok=True)
        
        # Remove burn-in
        burn_in = MCMC_CONFIG['n_burn']
        chains_burned = chains[:, burn_in:, :]
        log_post_burned = log_post[:, burn_in:]
        
        n_chains = chains.shape[0]
        
        # Gelman-Rubin diagnostics
        gr_stats = self.gelman_rubin_diagnostic(chains_burned)
        
        print("\n" + "="*80)
        print("CONVERGENCE DIAGNOSTICS (Gelman-Rubin statistic)")
        print("="*80)
        for i, param_name in enumerate(self.param_names):
            status = "✓ Converged" if gr_stats[i] < 1.1 else "✗ Not converged"
            print(f"{param_name:15s}: {gr_stats[i]:.4f}  {status}")
        
        # Trace plots
        fig, axes = plt.subplots(self.n_params, 1, figsize=(14, 3*self.n_params))
        if self.n_params == 1:
            axes = [axes]
        
        for i, param_name in enumerate(self.param_names):
            for chain_id in range(n_chains):
                axes[i].plot(chains_burned[chain_id, :, i], alpha=0.7, label=f'Chain {chain_id}')
            axes[i].set_ylabel(param_name, fontweight='bold')
            axes[i].legend()
            axes[i].grid(True, alpha=0.3)
        
        axes[-1].set_xlabel('Iteration (after burn-in)', fontweight='bold')
        plt.tight_layout()
        plt.savefig(f'{save_dir}/trace_plots.png', dpi=300, bbox_inches='tight')
        print(f"\n✓ Saved: trace_plots.png")
        plt.close()
        
        # Posterior distributions
        fig, axes = plt.subplots(1, self.n_params, figsize=(5*self.n_params, 4))
        if self.n_params == 1:
            axes = [axes]
        
        for i, param_name in enumerate(self.param_names):
            # Combine all chains
            all_samples = chains_burned[:, :, i].flatten()
            
            axes[i].hist(all_samples, bins=50, density=True, alpha=0.7, color='blue', edgecolor='black')
            axes[i].set_xlabel(param_name, fontweight='bold')
            axes[i].set_ylabel('Posterior density', fontweight='bold')
            axes[i].grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(f'{save_dir}/posterior_distributions.png', dpi=300, bbox_inches='tight')
        print(f"✓ Saved: posterior_distributions.png")
        plt.close()
        
        # Log-posterior evolution
        fig, ax = plt.subplots(figsize=(12, 5))
        for chain_id in range(n_chains):
            ax.plot(log_post_burned[chain_id], alpha=0.7, label=f'Chain {chain_id}')
        ax.set_xlabel('Iteration (after burn-in)', fontweight='bold')
        ax.set_ylabel('Log-posterior', fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(f'{save_dir}/log_posterior.png', dpi=300, bbox_inches='tight')
        print(f"✓ Saved: log_posterior.png")
        plt.close()
    
    def extract_posterior_estimates(self, chains):
        """Extract posterior mean and credible intervals"""
        burn_in = MCMC_CONFIG['n_burn']
        chains_burned = chains[:, burn_in:, :]
        
        # Combine all chains
        all_samples = chains_burned.reshape(-1, self.n_params)
        
        estimates = {}
        for i, param_name in enumerate(self.param_names):
            samples = all_samples[:, i]
            estimates[param_name] = {
                'mean': np.mean(samples),
                'median': np.median(samples),
                'std': np.std(samples),
                'ci_lower': np.percentile(samples, 2.5),
                'ci_upper': np.percentile(samples, 97.5),
            }
        
        return estimates
    
    def save_results(self, chains, estimates, save_dir='data/mcmc'):
        """Save estimated parameters to CSV"""
        os.makedirs(save_dir, exist_ok=True)
        
        # Create results dataframe
        results_df = pd.DataFrame({
            'Parameter': estimates.keys(),
            'Mean': [estimates[p]['mean'] for p in estimates.keys()],
            'Median': [estimates[p]['median'] for p in estimates.keys()],
            'Std': [estimates[p]['std'] for p in estimates.keys()],
            'CI_Lower_2.5%': [estimates[p]['ci_lower'] for p in estimates.keys()],
            'CI_Upper_97.5%': [estimates[p]['ci_upper'] for p in estimates.keys()],
        })
        
        results_df.to_csv(f'{save_dir}/mcmc_estimates.csv', index=False)
        print(f"✓ Saved: mcmc_estimates.csv")
        
        # Save full posterior samples
        burn_in = MCMC_CONFIG['n_burn']
        all_samples = chains[:, burn_in:, :].reshape(-1, self.n_params)
        
        samples_df = pd.DataFrame(all_samples, columns=self.param_names)
        samples_df.to_csv(f'{save_dir}/posterior_samples.csv', index=False)
        print(f"✓ Saved: posterior_samples.csv")
        
        return results_df

def main():
    # Setup
    unknown_params_dict = {name: idx for name, idx in UNKNOWN_PARAMS}
    likelihood = CancerCellDataLikelihood(EXPERIMENTAL_DATA, unknown_params_dict)
    sampler = MCMCDRAMSampler(likelihood, unknown_params_dict, PRIORS, PARAM_BOUNDS)
    
    # Run MCMC
    chains, log_post = sampler.run_parallel_chains()
    
    # Diagnostics and plots
    print("\n" + "="*80)
    print("CREATING DIAGNOSTIC PLOTS")
    print("="*80)
    sampler.plot_diagnostics(chains, log_post)
    
    # Extract and save results
    print("\n" + "="*80)
    print("SAVING RESULTS")
    print("="*80)
    estimates = sampler.extract_posterior_estimates(chains)
    results_df = sampler.save_results(chains, estimates)
    
    # Print summary
    print("\n" + "="*80)
    print("ESTIMATED PARAMETERS")
    print("="*80)
    print(results_df.to_string(index=False))
    print("\n" + "="*80)

if __name__ == '__main__':
    main()