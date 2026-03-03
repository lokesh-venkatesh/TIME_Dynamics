# TIMEdynamics: Bifurcation and Parameter Sensitivity Analysis Framework

A Julia-based framework for simulating, analyzing bifurcations, and estimating parameters in dynamical systems. Built with the 2D Lotka-Volterra model as a foundational example.

## Overview

This project provides a modular framework to:

1. **Simulate ODE systems** - Run trajectories for fixed time or until convergence
2. **Analyze bifurcations** - Systematically vary parameters and characterize attractors
3. **Estimate parameters** - Use Bayesian inference with MCMC to infer unknown parameters from data

## Project Structure

```
TIMEdynamics/
├── src/
│   ├── ode_system.jl          # ODE system definition (Lotka-Volterra)
│   └── parameters.jl          # Parameter sets and utilities
├── scripts/
│   ├── simulate.jl            # Forward simulation script
│   ├── bifurcation_analysis.jl # Bifurcation analysis (1D & 2D)
│   └── parameter_estimation.jl # Bayesian parameter estimation (MCMC)
├── data/                      # Output directory for results
└── README.md                  # This file
```

## Installation & Setup

### Prerequisites

- Julia 1.9+ (https://julialang.org/downloads/)
- Required packages (see below)

### Install Required Packages

Open Julia and run:

```julia
import Pkg
Pkg.add("DifferentialEquations")
Pkg.add("CSV")
Pkg.add("DataFrames")
Pkg.add("ArgParse")
Pkg.add("Statistics")
Pkg.add("LinearAlgebra")
Pkg.add("Distributions")
Pkg.add("StatsBase")
Pkg.add("Random")
Pkg.add("ProgressMeter")
```

Or use the provided `install_packages.jl` script:

```bash
julia install_packages.jl
```

## Quick Start

### 1. Run a Forward Simulation

```bash
# Simulate for a fixed time period
julia scripts/simulate.jl --time 100 --x0 1.0 0.5 --p default

# Simulate until the system settles to an attractor
julia scripts/simulate.jl --settle --x0 1.0 0.5 --p default

# Custom parameters
julia scripts/simulate.jl --time 50 --x0 2 1 --α 1.2 --β 0.12 --γ 0.4 --δ 0.015
```

**Output**: CSV file in `data/` directory containing time, prey, and predator columns.

### 2. Perform Bifurcation Analysis

#### 1D Bifurcation (vary one parameter)

```bash
# Vary alpha (prey growth rate) and detect how attractors change
julia scripts/bifurcation_analysis.jl --param α --range 0.5 2.0 --steps 50 --x0 1 0.5

# Vary beta (predation efficiency) with finer resolution
julia scripts/bifurcation_analysis.jl --param β --range 0.05 0.2 --steps 100

# Keep gamma fixed at a specific value
julia scripts/bifurcation_analysis.jl --param α --range 0.5 2.0 --steps 50 --fixed γ 0.4
```

#### 2D Bifurcation (vary two parameters)

```bash
# Create a 2D bifurcation diagram in (alpha, beta) space
julia scripts/bifurcation_analysis.jl \
    --param2d α β \
    --range1 0.5 2.0 --range2 0.05 0.2 \
    --steps1 40 --steps2 40

# Finer resolution
julia scripts/bifurcation_analysis.jl \
    --param2d α β \
    --range1 0.5 2.0 --range2 0.05 0.2 \
    --steps1 100 --steps2 100
```

**Output**: CSV file with parameter values and attractor characteristics.

### 3. Estimate Unknown Parameters

First, generate synthetic observation data:

```bash
julia scripts/simulate.jl --time 100 --x0 1.0 0.5 --α 1.3 --β 0.12 --γ 0.35 --δ 0.011 \
    -o synthetic_observations.csv
```

Then estimate the parameters using MCMC:

```bash
# Estimate alpha and beta
julia scripts/parameter_estimation.jl \
    --data data/synthetic_observations.csv \
    --param α β \
    --priors 0.5:2.0 0.05:0.2 \
    --iterations 5000 --burnin 1000 \
    --sigma 0.05

# Estimate multiple parameters with multiple chains
julia scripts/parameter_estimation.jl \
    --data data/synthetic_observations.csv \
    --param α β γ \
    --priors 0.5:2.0 0.05:0.2 0.1:0.5 \
    --iterations 10000 --burnin 2000 \
    --chains 4 --sigma 0.05
```

**Output**: 
- CSV file with posterior samples (MCMC chain)
- `.diag` file with diagnostics and summary statistics

## System Model: 2D Lotka-Volterra

The system implemented is the classic predator-prey (Lotka-Volterra) model:

$$\frac{dx}{dt} = \alpha x - \beta xy$$
$$\frac{dy}{dt} = -\gamma y + \delta xy$$

Where:
- **x** = prey population
- **y** = predator population
- **α** = intrinsic prey growth rate
- **β** = predation efficiency
- **γ** = predator mortality rate
- **δ** = conversion efficiency (new predators per prey consumed)

### Dynamics

- **Extinction equilibrium**: $(0, 0)$ - unstable
- **Coexistence equilibrium**: $\left(\frac{\gamma}{\delta}, \frac{\alpha}{\beta}\right)$ - center (neutrally stable)
- **Behavior**: The system exhibits **periodic oscillations** (limit cycles) for positive parameters

## File Format Specifications

### Simulation Output (CSV)

```
time,prey,predator
0.0,1.0,0.5
0.1,1.05,0.48
...
```

### Bifurcation Results (CSV)

**1D Analysis:**
```
param_value,attractor_type,x_mean,y_mean,x_std,y_std,dist_to_eq
0.5,limit_cycle,2.0,0.8,0.15,0.06,0.25
...
```

**2D Analysis:**
```
param1_value,param2_value,attractor_type,x_mean,y_mean,x_std,y_std,dist_to_eq
0.5,0.05,limit_cycle,2.0,0.8,0.15,0.06,0.25
...
```

### MCMC Output (CSV)

```
α,β,γ,δ
0.95,0.108,0.31,0.0099
0.97,0.111,0.32,0.0101
...
```

### Diagnostics File (.diag)

Text file containing acceptance rate, chain statistics, and posterior summary.

## Advanced Usage

### Custom ODE Systems

To modify the system or add new models:

1. Edit `src/ode_system.jl` - add your ODE function
2. Update `src/parameters.jl` - define parameter sets
3. Update scripts to call your custom ODE function

### Customizing Priors

For parameter estimation with non-uniform priors, modify the `parse_prior_specification()` function in `parameter_estimation.jl`.

### Parallel Execution

Run multiple bifurcation analyses or MCMC chains in parallel using Julia's built-in parallelism:

```bash
julia -p 4 scripts/bifurcation_analysis.jl ...
```

## Troubleshooting

**Simulation fails to converge**: Increase `--max-time` or reduce ODE solver tolerances in scripts.

**MCMC acceptance rate too low**: 
- Increase proposal scaling in `adaptive_metropolis_sampling()`
- Adjust priors to be more informative
- Check data for inconsistencies with the model

**Numerical instability**: 
- Reduce time step (increase saveat frequency)
- Tighten ODE solver tolerances (abstol, reltol)

## References

- **Lotka-Volterra Model**: Murray, J. D. (2003). Mathematical Biology I. Springer
- **MCMC Methods**: Haario et al. (2006). DRAM: A two-stage adaptive Metropolis algorithm. Statistics and Computing
- **Bifurcation Theory**: Kuznetsov, Y. A. (2004). Elements of Applied Bifurcation Theory. Springer

## Contributing

To extend this framework:

1. Add new ODE systems in `src/`
2. Implement specialized analysis methods in `scripts/`
3. Add test cases and validations

## License

Please specify your license here.

## Contact

For questions or issues, please contact the repository maintainers.
