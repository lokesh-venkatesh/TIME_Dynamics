# =============================================================================
# EXAMPLE 1: Basic Simulation
# =============================================================================

"""
import Pkg
Pkg.activate(".")


# Simulate for 100 time units with default parameters
julia> save("scripts/simulate.jl")
julia> run(`julia scripts/simulate.jl --time 100 --x0 1.0 0.5 --p default -v`)


# Simulate until system settles
julia> run(`julia scripts/simulate.jl --settle --x0 1.0 0.5 --p default -v`)

# Simulate with custom parameters
julia> run(`julia scripts/simulate.jl --time 50 --x0 2 1 --α 1.2 --β 0.12 --γ 0.4 --δ 0.015 -v`)

Output files will be created in: data/
"""