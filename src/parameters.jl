""" parameters.jl
Defines parameter sets and utilities for the Lotka-Volterra system. """

""" DEFAULT_PARAMETERS
Default parameter set for the 2D Lotka-Volterra system.
- α: intrinsic prey growth rate
- β: predation efficiency
- γ: predator mortality rate
- δ: conversion efficiency """
const DEFAULT_PARAMETERS = (
    α = 1.0,
    β = 0.1,
    γ = 0.3,
    δ = 0.01)

""" get_parameter_set(name::String)
Retrieve a pre-defined parameter set by name.
Arguments:  name::String: Name of the parameter set
Returns:    NamedTuple: Parameter set
Available sets: - "default": Default oscillatory behavior
                - "fast_predator": Fast predator dynamics
                - "strong_predation": Increased predation pressure """
function get_parameter_set(name::String)
    sets = Dict(
        "default" => DEFAULT_PARAMETERS,
        "fast_predator" => (α=1.0, β=0.1, γ=0.5, δ=0.01),
        "strong_predation" => (α=1.0, β=0.15, γ=0.3, δ=0.01),
        "weak_predation" => (α=1.0, β=0.05, γ=0.3, δ=0.01),
    )    
    if !haskey(sets, name)
        error("Unknown parameter set: $name. Available: $(join(keys(sets), ", "))")
    end
    return sets[name]
end

""" create_parameter_set(; α=1.0, β=0.1, γ=0.3, δ=0.01)
Create a custom parameter set.
Keyword Arguments:  α::Float64: Prey growth rate
                    β::Float64: Predation efficiency
                    γ::Float64: Predator mortality rate
                    δ::Float64: Conversion efficiency
Returns: NamedTuple: Custom parameter set """
function create_parameter_set(; α=1.0, β=0.1, γ=0.3, δ=0.01)
    return (α=α, β=β, γ=γ, δ=δ)
end

""" get_parameter_combinations(param_ranges::Dict)
Generate a grid of parameter combinations from specified ranges.
Arguments: param_ranges::Dict: Dictionary with parameter names as keys and ranges as values
                        e.g., Dict("α" => 0.5:0.1:1.5, "β" => 0.05:0.02:0.2)
Returns: Vector: Vector of parameter NamedTuples covering the grid
Example:    ranges = Dict("α" => 0.5:0.25:1.5, "β" => 0.05:0.05:0.15)
            combos = get_parameter_combinations(ranges)"""
function get_parameter_combinations(param_ranges::Dict)
    param_names = sort(collect(keys(param_ranges)))
    param_lists = [collect(param_ranges[name]) for name in param_names]
    combinations = vec(collect(Iterators.product(param_lists...))) # Create all combinations
    # Convert to NamedTuple format, replacing only specified parameters
    base_params = DEFAULT_PARAMETERS
    params_list = []
    
    for combo in combinations
        params_dict = Dict(param_names .=> combo)
        merged_params = merge(base_params, (; zip(Symbol.(param_names), combo)...))
        push!(params_list, merged_params)
    end
    return params_list
end

""" print_parameter_set(p)
Print a formatted parameter set.
Arguments: p::NamedTuple: Parameter set """
function print_parameter_set(p)
    println("Parameter Set:")
    for (key, val) in pairs(p)
        println("  $key = $val")
    end
end

export DEFAULT_PARAMETERS, get_parameter_set, create_parameter_set, 
       get_parameter_combinations, print_parameter_set