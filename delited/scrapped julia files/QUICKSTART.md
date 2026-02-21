# QUICKSTART Guide for TIMEdynamics

Get up and running in 5 minutes!

## 1. Install Dependencies (1 minute)

```bash
julia install_packages.jl
```

This installs all required Julia packages. You only need to do this once.

## 2. Run a Simple Simulation (1 minute)

```bash
julia scripts/simulate.jl --time 100 --x0 1.0 0.5 --p default --verbose
```

This simulates the Lotka-Volterra system for 100 time units starting from (x₀=1, y₀=0.5) with default parameters. A CSV file will be saved in `data/`.

## 3. Visualize the Trajectory (optional)

```bash
julia scripts/trajectory_analysis.jl --file data/trajectory_fixed_*.csv --verbose
```

This provides statistics about your trajectory.

## 4. Run 1D Bifurcation Analysis (2-3 minutes)

```bash
julia scripts/bifurcation_analysis.jl --param α --range 0.5 2.0 --steps 50 --verbose
```

This varies the prey growth rate (α) and identifies how the system's behavior changes. Results are saved as CSV.

## 5. Try Simulation Until Equilibrium

```bash
julia scripts/simulate.jl --settle --x0 1.0 0.5 --verbose
```

The system will run until it reaches an attractor (fixed point or limit cycle).

## 6. Try Parameter Estimation

First, generate synthetic data:
```bash
julia scripts/simulate.jl --time 100 --x0 1.0 0.5 --α 1.3 --β 0.12 --γ 0.35 --δ 0.011 -o synthetic_obs.csv
```

Then estimate parameters:
```bash
julia scripts/parameter_estimation.jl \
    --data data/synthetic_obs.csv \
    --param α β \
    --priors 0.5:2.0 0.05:0.2 \
    --iterations 5000 --burnin 1000 \
    --sigma 0.05 --verbose
```

## Common Commands

### Simulations

```bash
# Fixed time, default parameters
julia scripts/simulate.jl --time 100

# Until equilibrium, custom initial conditions
julia scripts/simulate.jl --settle --x0 2.0 1.0

# Custom parameters
julia scripts/simulate.jl --time 50 --α 1.2 --β 0.15 --γ 0.3 --δ 0.012

# Different parameter sets
julia scripts/simulate.jl --time 100 --p fast_predator
julia scripts/simulate.jl --time 100 --p strong_predation
```

### Bifurcation Analysis

```bash
# 1D: vary alpha
julia scripts/bifurcation_analysis.jl --param α --range 0.5 2.0 --steps 50

# 1D: vary beta with finer resolution
julia scripts/bifurcation_analysis.jl --param β --range 0.05 0.25 --steps 100

# 2D: vary alpha and beta together
julia scripts/bifurcation_analysis.jl --param2d α β --range1 0.5 2.0 --range2 0.05 0.2 --steps1 30 --steps2 30
```

### Parameter Estimation

```bash
# Estimate two parameters
julia scripts/parameter_estimation.jl --data data/synthetic_obs.csv --param α β --iterations 5000

# Estimate three parameters with multiple chains
julia scripts/parameter_estimation.jl --data data/synthetic_obs.csv --param α β γ \
    --iterations 10000 --burnin 2000 --chains 4
```

## Output Files

All results are saved in the `data/` directory:

- `trajectory_*.csv` - Simulation trajectories (time, prey, predator)
- `bifurcation_1d_*.csv` - 1D bifurcation results
- `bifurcation_2d_*.csv` - 2D bifurcation results
- `posterior_*.csv` - MCMC posterior samples (parameter estimation)
- `*.csv.diag` - Diagnostic information
- `*.csv.meta` - Metadata about the analysis

## Next Steps

1. **Read the full README.md** for detailed documentation
2. **Modify parameters** in `src/parameters.jl` to explore different behavior
3. **Implement your own ODE system** by editing `src/ode_system.jl`
4. **Analyze results** using your favorite tools (Python, MATLAB, Julia, etc.) - the CSV format is universal

## Help

For detailed options on any script:

```bash
julia scripts/simulate.jl --help
julia scripts/bifurcation_analysis.jl --help
julia scripts/parameter_estimation.jl --help
julia scripts/trajectory_analysis.jl --help
```

## Troubleshooting

**"Module not found" error?**
- Make sure you've run `julia install_packages.jl` first

**Simulation is slow?**
- Bifurcation analysis takes time as it runs many simulations
- Start with fewer steps (`--steps 20`) to test
- Use 2D analysis with coarser resolution first (`--steps1 20 --steps2 20`)

**MCMC acceptance rate very low?**
- Your priors might be too wide - make them more informative (closer to true values)
- Increase iterations and burn-in
- Check that your data makes sense for the model

**Questions?**
- See examples.jl for workflow demonstrations
- Check the full README.md for background theory and model details
