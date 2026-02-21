"""interactive_analysis.jl
Interactive Julia script for using the TIMEdynamics framework.
Run with: julia
Then at the REPL: include("interactive_analysis.jl")
This provides interactive functions for all major analysis types."""

push!(LOAD_PATH, joinpath(@__DIR__, "src"))
include(joinpath(@__DIR__, "src", "ode_system.jl"))
include(joinpath(@__DIR__, "src", "parameters.jl"))
using DifferentialEquations, CSV, DataFrames
using Statistics, LinearAlgebra, Distributions

# SIMULATION FUNCTIONS

""" sim_fixed(time_span::Float64, u0=nothing, param_set="default"; kwargs...)
Run a simulation for a fixed time period.
Examples:
    result = sim_fixed(100)                           # Default
    result = sim_fixed(100, [1.0, 0.5], "default")   # Custom IC
    result = sim_fixed(50, [2.0, 1.0], "fast_predator") """
function sim_fixed(time_span::Float64, u0=nothing, param_set="default"; kwargs...)
    if u0 === nothing
        u0 = [1.0, 0.5]
    end
    p = get_parameter_set(param_set)
    prob = ODEProblem(lotka_volterra!, u0, (0.0, time_span), p)
    sol = solve(prob, Tsit5(); saveat=0.1, kwargs...)
    return sol
end

""" sim_until_settle(u0=nothing, param_set="default"; max_time=1000.0)
Run a simulation until the system settles to an attractor.
Examples:
    result = sim_until_settle()
    result = sim_until_settle([2.0, 1.0], "fast_predator") """
function sim_until_settle(u0=nothing, param_set="default"; max_time=1000.0)
    if u0 === nothing
        u0 = [1.0, 0.5]
    end
    p = get_parameter_set(param_set)
    converged = [false] # Simple convergence callback
    history = Dict("states" => SVector[], "times" => Float64[])
    function condition(u, t, integrator)
        if length(history["states"]) < 10
            push!(history["states"], copy(u))
            push!(history["times"], t)
            return false
        else # Check if recent states are similar
            recent = history["states"][end-5:end]
            diffs = [norm(recent[i] - recent[1]) for i in 1:length(recent)]
            if all(d < 1e-3 for d in diffs)
                return true
            end
            if length(history["states"]) > 20 # Maintain recent history
                popfirst!(history["states"])
                popfirst!(history["times"])
            end
            push!(history["states"], copy(u))
            push!(history["times"], t)
            return false
        end
    end
    function affect!(integrator)
        terminate!(integrator)
    end
    cb = DiscreteCallback(condition, affect!, save_positions=(true, true))
    prob = ODEProblem(lotka_volterra!, u0, (0.0, max_time), p)
    sol = solve(prob, Tsit5(); callback=cb, saveat=0.1)
    return sol
end

# ANALYSIS FUNCTIONS

""" analyze_sol(sol)
Get statistics from a solution.
Returns: Dict with min, max, mean, std for prey and predator """
function analyze_sol(sol)
    x = sol[1, :]
    y = sol[2, :]
    return Dict( "x_min" => minimum(x), "x_max" => maximum(x),
        "x_mean" => mean(x), "x_std" => std(x), "y_min" => minimum(y),
        "y_max" => maximum(y), "y_mean" => mean(y), "y_std" => std(y))
end
""" attractor_type(sol)
Determine if trajectory is at fixed point or limit cycle.
Returns: "fixed_point" or "limit_cycle" """
function attractor_type(sol)
    x = sol[1, :]
    y = sol[2, :]
    if std(x) < 0.01 * mean(x) && std(y) < 0.01 * mean(y)
        return "fixed_point"
    else
        return "limit_cycle"
    end
end

""" bifurc_1d(param_name::String, param_range::Tuple, n_steps::Int; 
              u0=nothing, fixed_params=Dict())
Perform 1D bifurcation analysis.
Example:
    results = bifurc_1d("α", (0.5, 2.0), 50) """
function bifurc_1d(param_name::String, param_range::Tuple, n_steps::Int;
                   u0=nothing, fixed_params=Dict())
    if u0 === nothing
        u0 = [1.0, 0.5]
    end
    param_sym = Symbol(param_name)
    param_values = range(param_range[1], param_range[2], length=n_steps)
    results = DataFrame(param_value = Float64[], attractor_type = String[],
        x_mean = Float64[], y_mean = Float64[],
        x_std = Float64[], y_std = Float64[])
    p_base = merge(DEFAULT_PARAMETERS, fixed_params)
    print("Progress: ")
    for (i, p_val) in enumerate(param_values)
        print(".")
        p = merge(p_base, NamedTuple{(param_sym,)}((p_val,)))
        try
            prob = ODEProblem(lotka_volterra!, u0, (0.0, 500.0), p)
            sol = solve(prob, Tsit5(); saveat=0.1)
            stats = analyze_sol(sol) 
            push!(results, (p_val,attractor_type(sol),
                stats["x_mean"], stats["y_mean"],
                stats["x_std"], stats["y_std"]))
        catch e
            # Skip failed points
        end
    end
    println()
    return results
end

""" plot_trajectory(sol; title="Trajectory")
Plot the trajectory (requires Plots.jl) """
function plot_trajectory(sol; title="Trajectory")
    try
        using Plots
        t = sol.t
        x = sol[1, :]
        y = sol[2, :]
        p1 = plot(t, x, label="Prey", xlabel="Time")
        p2 = plot(t, y, label="Predator", xlabel="Time")
        p3 = plot(x, y, label="Phase portrait", xlabel="Prey", ylabel="Predator")
        return plot(p1, p2, p3, layout=(3,1), title=title)
    catch
        println("Plots.jl not available. Install with: Pkg.add(\"Plots\")")
    end
end

""" plot_bifurcation(results; param_name="α")
Plot 1D bifurcation diagram. """
function plot_bifurcation(results; param_name="α")
    try
        using Plots
        fixed_pts = results[results.attractor_type .== "fixed_point", :]
        limit_cyc = results[results.attractor_type .== "limit_cycle", :]
        p1 = scatter(fixed_pts.param_value, fixed_pts.x_mean, label="Fixed point")
        scatter!(p1, limit_cyc.param_value, limit_cyc.x_mean .- limit_cyc.x_std, 
                 label="Limit cycle", alpha=0.5)
        scatter!(p1, limit_cyc.param_value, limit_cyc.x_mean .+ limit_cyc.x_std,
                 label="", alpha=0.5)
        xlabel!(p1, "$param_name")
        ylabel!(p1, "Prey population")
        return p1
    catch
        println("Plots.jl not available.")
    end
end

# WORKFLOWS

""" workflow_explore()
Interactive workflow: Explore system dynamics."""
function workflow_explore()
    println("""EXPLORATION WORKFLOW""")
    println("1. Run simulation with default parameters...")
    sol1 = sim_fixed(100)
    stats1 = analyze_sol(sol1)
    println("   Attractor type: $(attractor_type(sol1))")
    println("   Prey: $(round(stats1["x_mean"], digits=3)) ± $(round(stats1["x_std"], digits=3))")
    println("   Predator: $(round(stats1["y_mean"], digits=3)) ± $(round(stats1["y_std"], digits=3))")
    println()
    println("2. Try fast_predator set...")
    sol2 = sim_fixed(100, [1.0, 0.5], "fast_predator")
    stats2 = analyze_sol(sol2)
    println("   Attractor type: $(attractor_type(sol2))")
    println("   Prey: $(round(stats2["x_mean"], digits=3)) ± $(round(stats2["x_std"], digits=3))")
    println()
    println("3. Quick bifurcation check (vary α)...")
    bifurc = bifurc_1d("α", (0.5, 1.5), 20)
    println("   Transitions found:")
    for type in unique(bifurc.attractor_type)
        idxs = bifurc.attractor_type .== type
        param_range = (minimum(bifurc[idxs, 1]), maximum(bifurc[idxs, 1]))
        println("   - $type in α ∈ $(param_range)")
    end
end

"""
    workflow_bifurcation_survey()

Interactive workflow: Survey bifurcations.
"""
function workflow_bifurcation_survey()
    println("""BIFURCATION SURVEY WORKFLOW""")
    println("Performing 1D bifurcation analysis in α...")
    results = bifurc_1d("α", (0.3, 3.0), 100)
    println("\nFound $(length(unique(results.attractor_type))) attractor type(s)")
    println("Bifurcation diagram computed: $(nrow(results)) points")
    transitions = Int[] # Print transition points
    for i in 2:nrow(results)
        if results[i, 2] != results[i-1, 2]
            push!(transitions, i)
        end
    end
    if !isempty(transitions)
        println("\nBifurcation transitions occur near:")
        for i in transitions
            α_val = results[i, 1]
            println("  α ≈ $α_val")
        end
    end
    return results
end

# HELP & INFO

""" list_functions()
List all available interactive functions. """
function list_functions()
    println("""TIMEdynamics Interactive Functions
    
    SIMULATION
    ──────────
    • sim_fixed(time, u0, param_set)       - Fixed time simulation
    • sim_until_settle(u0, param_set)      - Simulate to attractor
    
    ANALYSIS
    ────────
    • analyze_sol(sol)                     - Get statistics
    • attractor_type(sol)                  - Classify attractor
    • bifurc_1d(param, range, steps)       - 1D bifurcation diagram
    
    VISUALIZATION (requires Plots.jl)
    ─────────────────────────────────
    • plot_trajectory(sol)                 - Plot trajectory
    • plot_bifurcation(results)            - Plot bifurcation diagram
    
    WORKFLOWS
    ─────────
    • workflow_explore()                   - Explore dynamics
    • workflow_bifurcation_survey()        - Survey bifurcations
    
    UTILITIES
    ─────────
    • list_functions()                     - Show this list
    • export_results(results, filename)    - Save to CSV
    
    EXAMPLES
    ────────
    julia> sol = sim_fixed(100, [1.0, 0.5], "default")
    julia> stats = analyze_sol(sol)
    julia> plot_trajectory(sol)
    
    julia> results = bifurc_1d("α", (0.5, 2.0), 50)
    julia> plot_bifurcation(results)
    
    julia> workflow_explore()
    """)
end

""" export_results(results::DataFrame, filename::String)
Save DataFrame to CSV. """
function export_results(results::DataFrame, filename::String)
    data_dir = joinpath(@__DIR__, "data")
    mkpath(data_dir)
    filepath = joinpath(data_dir, filename)
    CSV.write(filepath, results)
    println("Saved to: $filepath")
end

# STARTUP MESSAGE

println("""
╔═══════════════════════════════════════════════════════════════╗
║  TIMEdynamics - Interactive Analysis Session                 ║
╚═══════════════════════════════════════════════════════════════╝

Type: list_functions()  to see all available commands

Quick examples:

  sol = sim_fixed(100)
  stats = analyze_sol(sol)
  attractor_type(sol)
  
  results = bifurc_1d("α", (0.5, 2.0), 50)
  
  workflow_explore()

""")
