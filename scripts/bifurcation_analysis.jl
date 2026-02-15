""" bifurcation_analysis.jl
    Bifurcation analysis script for the Lotka-Volterra system.
Usage:  julia bifurcation_analysis.jl --param α --range 0.5 2.0 --steps 50
        julia bifurcation_analysis.jl --param β --range 0.05 0.2 --steps 100 --fixed α 1.2
        julia bifurcation_analysis.jl --param2d α β --range1 0.5 2.0 --range2 0.05 0.2 --steps1 30 --steps2 30 """

using DifferentialEquations
using CSV, DataFrames
using ArgParse
using Statistics
using LinearAlgebra
using Dates
using ProgressMeter

include(joinpath(@__DIR__, "..", "src", "ode_system.jl")) # Load modules directly
include(joinpath(@__DIR__, "..", "src", "parameters.jl"))

""" parse_commandline()
Parse command-line arguments for bifurcation analysis.
Returns: Dict: Parsed arguments """
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--param", "-p"
            help = "parameter to vary (α, β, γ, δ)"
            arg_type = String
        "--param2d"
            help = "two parameters to vary for 2D bifurcation diagram"
            nargs = 2
            arg_type = String
        "--range"
            help = "parameter range [min max] for 1D analysis"
            nargs = 2
            arg_type = Float64
        "--range1"
            help = "first parameter range [min max] for 2D analysis"
            nargs = 2
            arg_type = Float64
        "--range2"
            help = "second parameter range [min max] for 2D analysis"
            nargs = 2
            arg_type = Float64
        "--steps", "-n"
            help = "number of parameter values to sample"
            arg_type = Int
            default = 50
        "--steps1"
            help = "number of samples for first parameter (2D)"
            arg_type = Int
            default = 30
        "--steps2"
            help = "number of samples for second parameter (2D)"
            arg_type = Int
            default = 30
        "--fixed"
            help = "fixed parameter values [param_name param_value ...]"
            nargs = '+'
            arg_type = String
        "--x0"
            help = "initial conditions [x_prey x_predator]"
            nargs = 2
            arg_type = Float64
            default = [1.0, 0.5]
        "--output", "-o"
            help = "output filename (in data/ directory)"
            arg_type = String
            default = nothing
        "--max-time"
            help = "maximum simulation time"
            arg_type = Float64
            default = 1000.0
        "--transient"
            help = "time to wait before analysis (let transients die out)"
            arg_type = Float64
            default = 100.0
        "--verbose", "-v"
            help = "print detailed output"
            action = :store_true
    end
    return parse_args(s)
end

""" parse_fixed_parameters(fixed_args::Vector{String})
Parse fixed parameter specifications from command line.
Arguments: fixed_args::Vector{String}: Vector of ["param_name", "param_value", ...]
Returns: Dict: Dictionary of parameter name -> value """
function parse_fixed_parameters(fixed_args::Vector{String})
    fixed_params = Dict()
    for i in 1:2:length(fixed_args)
        param_name = Symbol(fixed_args[i])
        param_value = parse(Float64, fixed_args[i+1])
        fixed_params[param_name] = param_value
    end
    return fixed_params
end

""" detect_attractor(sol, p; transient_idx=1)
Detect and characterize the attractor from a solution.
Arguments:  sol::ODESolution: Trajectory solution
            p::NamedTuple: Parameters
            transient_idx::Int: Index after which to analyze (skip transients)
Returns: Dict: Attractor properties (type, equilibrium, period, etc.)
"""
function detect_attractor(sol, p; transient_idx=1)
    result = Dict()
    t = sol.t[transient_idx:end] # Skip transient part
    x = sol[1, transient_idx:end]
    y = sol[2, transient_idx:end]
    x_mean = mean(x) # Calculate amplitude and mean
    y_mean = mean(y)
    x_std = std(x)
    y_std = std(y)
    result["x_mean"] = x_mean
    result["y_mean"] = y_mean
    result["x_std"] = x_std
    result["y_std"] = y_std
    result["x_min"] = minimum(x)
    result["x_max"] = maximum(x)
    result["y_min"] = minimum(y)
    result["y_max"] = maximum(y)
    if x_std < 1e-4 && y_std < 1e-4 # Detect attractor type
        result["type"] = "fixed_point"
    else
        result["type"] = "limit_cycle"
        peaks_x = findpeaks(x) # Try to estimate period by finding peaks
        if length(peaks_x) > 1
            time_diff = diff(t[peaks_x])
            if length(time_diff) > 0
                result["period_estimate"] = mean(time_diff)
            end
        end
    end
    eq_unstable, eq_center = get_equilibrium_points(p) # Compare to theoretical equilibrium
    dist_to_center = sqrt((x_mean - eq_center[1])^2 + (y_mean - eq_center[2])^2)
    result["dist_to_center_eq"] = dist_to_center
    return result
end

""" findpeaks(x; threshold=0.1)
Simple peak detection in a time series.
Arguments:  x::Vector: Time series data
            threshold::Float64: Relative threshold for peak detection
Returns:    Vector{Int}: Indices of detected peaks """
function findpeaks(x; threshold=0.1)
    peaks = Int[]
    for i in 2:length(x)-1
        if x[i] > x[i-1] && x[i] > x[i+1]
            push!(peaks, i)
        end
    end
    return peaks
end

""" run_bifurcation_1d(param_name, param_range, n_steps, u0, fixed_params)
Perform 1D bifurcation analysis along a single parameter.
Arguments:  param_name::String: Parameter to vary
            param_range::Tuple: (min, max) of parameter range
            n_steps::Int: Number of parameter values
            u0::Vector: Initial conditions
            fixed_params::Dict: Other fixed parameters
Returns: DataFrame: Results with columns [param, attractor_type, x_mean, y_mean, period, ...] """
function run_bifurcation_1d(param_name, param_range, n_steps, u0, fixed_params, max_time, transient_time)
    param_sym = Symbol(param_name)
    param_values = range(param_range[1], param_range[2], length=n_steps)
    results = DataFrame(param_value = Float64[], attractor_type = String[],
        x_mean = Float64[], y_mean = Float64[], x_std = Float64[],
        y_std = Float64[], x_min = Float64[], x_max = Float64[],
        y_min = Float64[], y_max = Float64[], dist_to_eq = Float64[])
    @showprogress for p_val in param_values
        p = merge(DEFAULT_PARAMETERS, fixed_params) # Create parameter set
        p = merge(p, Dict(param_sym => p_val))
        try # Run simulation
            prob = ODEProblem(lotka_volterra!, u0, (0.0, max_time), p)
            sol = solve(prob, Tsit5(), saveat=0.1, abstol=1e-8, reltol=1e-8)
            transient_idx = max(1, Int(floor(transient_time / 0.1))) # Detect attractor
            attractor = detect_attractor(sol, p; transient_idx=transient_idx)
            push!(results, (p_val, # Store results
                attractor["type"], attractor["x_mean"], attractor["y_mean"],
                attractor["x_std"], attractor["y_std"], attractor["x_min"],
                attractor["x_max"], attractor["y_min"], attractor["y_max"],
                attractor["dist_to_center_eq"]))
        catch e
            if verbose
                println("Error at $param_name = $p_val: $e")
            end
        end
    end
    return results
end

""" run_bifurcation_2d(param1_name, param2_name, param1_range, param2_range, 
                       n_steps1, n_steps2, u0, fixed_params)
Perform 2D bifurcation analysis over two parameters.
Arguments:  param1_name::String: First parameter to vary
            param2_name::String: Second parameter to vary
            param1_range::Tuple: (min, max) of first parameter range
            param2_range::Tuple: (min, max) of second parameter range
            n_steps1::Int: Number of values for first parameter
            n_steps2::Int: Number of values for second parameter
            u0::Vector: Initial conditions
            fixed_params::Dict: Other fixed parameters
Returns: DataFrame: Results with columns [param1, param2, attractor_type, x_mean, y_mean, ...] """
function run_bifurcation_2d(param1_name, param2_name, param1_range, param2_range,
                            n_steps1, n_steps2, u0, fixed_params, max_time, transient_time)
    param1_sym = Symbol(param1_name)
    param2_sym = Symbol(param2_name)
    param1_values = range(param1_range[1], param1_range[2], length=n_steps1)
    param2_values = range(param2_range[1], param2_range[2], length=n_steps2)
    results = DataFrame(param1_value = Float64[], param2_value = Float64[],
        attractor_type = String[], x_mean = Float64[], y_mean = Float64[],
        x_std = Float64[], y_std = Float64[], dist_to_eq = Float64[])
    
    n_total = n_steps1 * n_steps2
    @showprogress dt=1.0 for p1_val in param1_values, p2_val in param2_values
        p = merge(DEFAULT_PARAMETERS, fixed_params) # Create parameter set
        p = merge(p, Dict(param1_sym => p1_val, param2_sym => p2_val))
        try # Run simulation
            prob = ODEProblem(lotka_volterra!, u0, (0.0, max_time), p)
            sol = solve(prob, Tsit5(), saveat=0.1, abstol=1e-8, reltol=1e-8)
            transient_idx = max(1, Int(floor(transient_time / 0.1))) # Detect attractor
            attractor = detect_attractor(sol, p; transient_idx=transient_idx)
            push!(results, (p1_val, p2_val, # Store results
                attractor["type"], attractor["x_mean"], attractor["y_mean"],
                attractor["x_std"], attractor["y_std"], attractor["dist_to_center_eq"]))
        catch e # Skip points that fail
        end
    end
    return results
end

""" save_bifurcation_results(results, param_names, filename)
Save bifurcation analysis results to file.
Arguments:  results::DataFrame: Bifurcation results
            param_names::Vector{String}: Names of varied parameters
            filename::String: Output filename (without directory)"""
function save_bifurcation_results(results, param_names, filename)
    data_dir = joinpath(@__DIR__, "..", "data")
    mkpath(data_dir)
    filepath = joinpath(data_dir, filename)
    CSV.write(filepath, results)
    metadata_file = filepath * ".meta" # Save metadata
    open(metadata_file, "w") do f
        write(f, "Bifurcation Analysis Metadata\n")
        write(f, "=============================\n")
        write(f, "Generated: $(now())\n")
        write(f, "Varied parameters: $(join(param_names, ", "))\n")
        write(f, "Number of points analyzed: $(nrow(results))\n")
    end
    return filepath
end

""" main()
    Main execution function. """
function main()
    args = parse_commandline()
    if args["param"] === nothing && args["param2d"] === nothing # Validate arguments
        error("Must specify either --param (1D) or --param2d (2D) for bifurcation analysis")
    end
    fixed_params_dict = Dict() # Parse fixed parameters
    if args["fixed"] !== nothing
        fixed_parsed = parse_fixed_parameters(args["fixed"])
        fixed_params_dict = fixed_parsed
    end
    u0 = args["x0"]
    max_time = args["max-time"]
    transient_time = args["transient"]
    println("=" ^ 70)
    println("BIFURCATION ANALYSIS")
    println("=" ^ 70)
    println("Initial conditions: x₀ = $(u0[1]), y₀ = $(u0[2])")
    println()
    if args["param"] !== nothing # Run 1D bifurcation
        param_name = args["param"]
        param_range = tuple(args["range"]...)
        n_steps = args["steps"]
        println("1D Bifurcation Analysis")
        println("Parameter: $param_name, Range: $(param_range[1]) to $(param_range[2])")
        println("Steps: $n_steps")
        println()
        results = run_bifurcation_1d(param_name, param_range, n_steps, u0, fixed_params_dict, max_time, transient_time)
        if args["output"] === nothing # Determine output filename
            timestamp = Dates.format(now(), "yyyymmdd_HHMMSS")
            output_filename = "bifurcation_1d_$(param_name)_$(timestamp).csv"
        else
            output_filename = args["output"]
        end
        filepath = save_bifurcation_results(results, [param_name], output_filename)
        println("Results saved to: $filepath")
        println("\nSummary:")
        println("  Attractor types found: $(unique(results.attractor_type))")
    elseif args["param2d"] !== nothing # Run 2D bifurcation
        param1_name = args["param2d"][1]
        param2_name = args["param2d"][2]
        param1_range = tuple(args["range1"]...)
        param2_range = tuple(args["range2"]...)
        n_steps1 = args["steps1"]
        n_steps2 = args["steps2"]
        println("2D Bifurcation Analysis")
        println("Parameters: $param1_name ($(param1_range[1]) to $(param1_range[2])), $param2_name ($(param2_range[1]) to $(param2_range[2]))")
        println("Steps: $n_steps1 × $n_steps2")
        println()
        results = run_bifurcation_2d(param1_name, param2_name, param1_range, param2_range,
                                     n_steps1, n_steps2, u0, fixed_params_dict, max_time, transient_time)
        if args["output"] === nothing # Determine output filename
            timestamp = Dates.format(now(), "yyyymmdd_HHMMSS")
            output_filename = "bifurcation_2d_$(param1_name)_$(param2_name)_$(timestamp).csv"
        else
            output_filename = args["output"]
        end
        filepath = save_bifurcation_results(results, [param1_name, param2_name], output_filename)
        println("Results saved to: $filepath")
        println("\nSummary:")
        println("  Attractor types found: $(unique(results.attractor_type))")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
