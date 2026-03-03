"""
    examples.jl

Example workflows demonstrating use of the TIMEdynamics framework.

This file contains Julia commands that can be run to demonstrate:
1. Running simulations
2. Performing bifurcation analysis
3. Parameter estimation

Run individual lines or sections at the Julia REPL or copy-paste as needed.
"""

import Pkg
Pkg.activate(".")

# =============================================================================
# EXAMPLE 1: Basic Simulation
# =============================================================================

println("""
╔════════════════════════════════════════════════════════════════╗
║  EXAMPLE 1: BASIC SIMULATIONS                                ║
╚════════════════════════════════════════════════════════════════╝

Try these commands:

# Simulate for 100 time units with default parameters
julia> save("scripts/simulate.jl")
julia> run(`julia scripts/simulate.jl --time 100 --x0 1.0 0.5 --p default -v`)

# Simulate until system settles
julia> run(`julia scripts/simulate.jl --settle --x0 1.0 0.5 --p default -v`)

# Simulate with custom parameters
julia> run(`julia scripts/simulate.jl --time 50 --x0 2 1 --α 1.2 --β 0.12 --γ 0.4 --δ 0.015 -v`)

Output files will be created in: data/

""")

# =============================================================================
# EXAMPLE 2: 1D Bifurcation Analysis
# =============================================================================

println("""
╔════════════════════════════════════════════════════════════════╗
║  EXAMPLE 2: 1D BIFURCATION ANALYSIS                          ║
╚════════════════════════════════════════════════════════════════╝

Try these commands:

# 1D bifurcation: vary alpha (prey growth rate)
julia> run(`julia scripts/bifurcation_analysis.jl --param α --range 0.5 2.0 --steps 50 --x0 1 0.5 -v`)

# 1D bifurcation: vary beta (predation efficiency) with finer resolution
julia> run(`julia scripts/bifurcation_analysis.jl --param β --range 0.05 0.25 --steps 100 -v`)

# 1D bifurcation: vary gamma with other parameters fixed
julia> run(`julia scripts/bifurcation_analysis.jl --param γ --range 0.1 0.5 --steps 50 --fixed α 1.0 β 0.1 -v`)

Output files: bifurcation_1d_*.csv

""")

# =============================================================================
# EXAMPLE 3: 2D Bifurcation Analysis
# =============================================================================

println("""
╔════════════════════════════════════════════════════════════════╗
║  EXAMPLE 3: 2D BIFURCATION ANALYSIS                          ║
╚════════════════════════════════════════════════════════════════╝

Try these commands:

# 2D bifurcation: vary alpha and beta together
# (This may take a few minutes)
julia> run(`julia scripts/bifurcation_analysis.jl --param2d α β --range1 0.5 2.0 --range2 0.05 0.25 --steps1 30 --steps2 30 -v`)

# Higher resolution (will take longer)
julia> run(`julia scripts/bifurcation_analysis.jl --param2d α β --range1 0.5 2.0 --range2 0.05 0.25 --steps1 50 --steps2 50 -v`)

Output files: bifurcation_2d_*.csv

""")

# =============================================================================
# EXAMPLE 4: Parameter Estimation
# =============================================================================

println("""
╔════════════════════════════════════════════════════════════════╗
║  EXAMPLE 4: PARAMETER ESTIMATION                             ║
╚════════════════════════════════════════════════════════════════╝

Step 1: Generate synthetic observation data
julia> run(`julia scripts/simulate.jl --time 100 --x0 1.0 0.5 --α 1.3 --β 0.12 --γ 0.35 --δ 0.011 -o synthetic_obs.csv`)

Step 2: Estimate the parameters using MCMC
julia> run(`julia scripts/parameter_estimation.jl \\
    --data data/synthetic_obs.csv \\
    --param α β \\
    --priors 0.5:2.0 0.05:0.2 \\
    --iterations 5000 --burnin 1000 \\
    --sigma 0.05 -v`)

Step 3: Estimate more parameters (requires more data)
julia> run(`julia scripts/parameter_estimation.jl \\
    --data data/synthetic_obs.csv \\
    --param α β γ δ \\
    --priors 0.5:2.0 0.05:0.2 0.1:0.5 0.005:0.02 \\
    --iterations 10000 --burnin 2000 \\
    --chains 4 --sigma 0.05 -v`)

Output files: posterior_*.csv
Diagnostics: posterior_*.csv.diag

""")

# =============================================================================
# EXAMPLE 5: Complete Workflow
# =============================================================================

println("""
╔════════════════════════════════════════════════════════════════╗
║  EXAMPLE 5: COMPLETE WORKFLOW                                  ║
╚════════════════════════════════════════════════════════════════╝

A typical analysis workflow:

1. EXPLORE: Run simulations with different parameters
   julia> run(`julia scripts/simulate.jl --time 200 --p default -v`)
   julia> run(`julia scripts/simulate.jl --time 200 --p fast_predator -v`)

2. IDENTIFY BIFURCATIONS: Perform 1D analysis
   julia> run(`julia scripts/bifurcation_analysis.jl --param α --range 0.5 2.0 --steps 50 -v`)

3. MAP PARAMETER SPACE: Perform 2D analysis
   julia> run(`julia scripts/bifurcation_analysis.jl --param2d α β \\
      --range1 0.5 2.0 --range2 0.05 0.25 --steps1 30 --steps2 30 -v`)

4. INFER PARAMETERS: Estimate from data if you have measurements
   (See Example 4)

5. ANALYZE RESULTS: Use your favorite plotting/analysis tools
   - Load the CSV files
   - Create bifurcation diagrams
   - Plot MCMC traces and posteriors

""")

# =============================================================================
# EXAMPLE 6: Analyzing Results
# =============================================================================

println("""
╔════════════════════════════════════════════════════════════════╗
║  EXAMPLE 6: ANALYZING RESULTS IN JULIA                       ║
╚════════════════════════════════════════════════════════════════╝

After running simulations, you can analyze results in Julia:

# Load a simulation
julia> using CSV, DataFrames
julia> traj = CSV.read("data/trajectory_fixed_*.csv", DataFrame)
julia> plot(traj.time, traj.prey, label="Prey")
julia> plot!(traj.time, traj.predator, label="Predator")

# Load bifurcation results
julia> bifurc = CSV.read("data/bifurcation_1d_α_*.csv", DataFrame)
julia> plot(bifurc.param_value, bifurc.x_mean, seriestype=:scatter)

# Load MCMC results
julia> posterior = CSV.read("data/posterior_α_β_*.csv", DataFrame)
julia> using StatsPlots
julia> @df posterior plot(:α, :β, seriestype=:scatter, alpha=0.3)

""")

println("""
═══════════════════════════════════════════════════════════════════

For full documentation, see README.md
For detailed parameter help on each script, use:
  julia scripts/simulate.jl --help
  julia scripts/bifurcation_analysis.jl --help
  julia scripts/parameter_estimation.jl --help

═══════════════════════════════════════════════════════════════════
""")
