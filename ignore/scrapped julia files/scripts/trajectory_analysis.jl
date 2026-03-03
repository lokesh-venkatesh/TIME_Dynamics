""" trajectory_analysis.jl
Utility script for analyzing simulation trajectories.
Usage:
    julia trajectory_analysis.jl --file trajectory_file.csv --output analysis_report.txt
    julia trajectory_analysis.jl --file trajectory_file.csv --plot --figure out.png"""

using CSV, DataFrames
using ArgParse
using Statistics
using Dates

push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
include(joinpath(@__DIR__, "..", "src", "ode_system.jl"))

""" parse_commandline()
Parse command-line arguments.
Returns: Dict: Parsed arguments """
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--file", "-f"
            help = "trajectory CSV file to analyze"
            arg_type = String
            required = true
        "--output", "-o"
            help = "output file for analysis report"
            arg_type = String
            default = nothing
        "--plot"
            help = "generate plots"
            action = :store_true
        "--figure"
            help = "save plots to file"
            arg_type = String
            default = nothing
        "--verbose", "-v"
            help = "verbose output"
            action = :store_true
    end
    return parse_args(s)
end

""" analyze_trajectory(df)
Perform comprehensive analysis of a trajectory.
Arguments: df::DataFrame: Trajectory data (columns: time, prey, predator)
Returns: Dict: Analysis results """
function analyze_trajectory(df)
    results = Dict()
    t = df.time
    x = df.prey
    y = df.predator
    results["n_points"] = length(t) # Basic statistics
    results["t_start"] = t[1]
    results["t_end"] = t[end]
    results["duration"] = t[end] - t[1] 
    results["x_min"] = minimum(x) # Prey statistics
    results["x_max"] = maximum(x)
    results["x_mean"] = mean(x)
    results["x_median"] = median(x)
    results["x_std"] = std(x)
    results["x_cv"] = results["x_std"] / results["x_mean"]  # Coefficient of variation
    results["y_min"] = minimum(y) # Predator statistics
    results["y_max"] = maximum(y)
    results["y_mean"] = mean(y)
    results["y_median"] = median(y)
    results["y_std"] = std(y)
    results["y_cv"] = results["y_std"] / results["y_mean"]
    results["x_final"] = x[end] # Final state
    results["y_final"] = y[end]
    # Detect if settled to fixed point or limit cycle
    if results["x_std"] < 0.01 * results["x_mean"] && results["y_std"] < 0.01 * results["y_mean"]
        results["attractor_type"] = "fixed_point"
    else
        results["attractor_type"] = "limit_cycle"
    end
    mid_idx = Int(floor(length(t) / 2)) # Look for transient behavior by splitting analysis
    x_first = x[1:mid_idx]
    x_second = x[mid_idx+1:end]
    y_first = y[1:mid_idx]
    y_second = y[mid_idx+1:end]
    results["transient_period"] = t[mid_idx]
    results["x_change_ratio"] = (mean(x_second) - mean(x_first)) / mean(x_first)
    results["y_change_ratio"] = (mean(y_second) - mean(y_first)) / mean(y_first)
    return results
end

""" print_analysis_report(results, df)
Print a formatted analysis report.
Arguments: results::Dict: Analysis results
            df::DataFrame: Trajectory data """
function print_analysis_report(results, df)
    report = """
    ═══════════════════════════════════════════════════════════════
    TRAJECTORY ANALYSIS REPORT
    ═══════════════════════════════════════════════════════════════
    
    SIMULATION DURATION
    ───────────────────
    Time range: $(results["t_start"]) to $(results["t_end"])
    Duration: $(results["duration"])
    Number of time points: $(results["n_points"])
    
    ATTRACTOR CHARACTERISTICS
    ─────────────────────────
    Type: $(results["attractor_type"])
    
    PREY (x) STATISTICS
    ───────────────────
    Minimum: $(round(results["x_min"], digits=6))
    Maximum: $(round(results["x_max"], digits=6))
    Mean: $(round(results["x_mean"], digits=6))
    Median: $(round(results["x_median"], digits=6))
    Std Dev: $(round(results["x_std"], digits=6))
    Coeff. of Variation: $(round(results["x_cv"], digits=4))
    Final Value: $(round(results["x_final"], digits=6))
    
    PREDATOR (y) STATISTICS
    ───────────────────────
    Minimum: $(round(results["y_min"], digits=6))
    Maximum: $(round(results["y_max"], digits=6))
    Mean: $(round(results["y_mean"], digits=6))
    Median: $(round(results["y_median"], digits=6))
    Std Dev: $(round(results["y_std"], digits=6))
    Coeff. of Variation: $(round(results["y_cv"], digits=4))
    Final Value: $(round(results["y_final"], digits=6))
    
    TRANSIENT ANALYSIS
    ──────────────────
    Transient Period: $(results["transient_period"])
    Prey Change (first half to second half): $(round(results["x_change_ratio"], digits=4))
    Predator Change (first half to second half): $(round(results["y_change_ratio"], digits=4))
    
    ═══════════════════════════════════════════════════════════════
    """
    return report
end

""" main()
    Main execution function. """
function main()
    args = parse_commandline()
    if !isfile(args["file"]) # Check if file exists
        error("File not found: $(args["file"])")
    end
    df = CSV.read(args["file"], DataFrame) # Load data
    results = analyze_trajectory(df) # Analyze
    report = print_analysis_report(results, df) # Print report
    println(report)
    if args["output"] !== nothing # Save report if requested
        open(args["output"], "w") do f
            write(f, report)
        end
        println("\nReport saved to: $(args["output"])")
    end
    if args["plot"] # Plot if requested
        try
            using Plots
            p1 = plot(df.time, df.prey, label="Prey (x)", xlabel="Time", ylabel="Population")
            p2 = plot(df.time, df.predator, label="Predator (y)", xlabel="Time", ylabel="Population")
            p3 = plot(df.prey, df.predator, label="Phase portrait", xlabel="Prey", ylabel="Predator", seriestype=:path)
            p = plot(p1, p2, p3, layout=(3,1), size=(800, 900))
            if args["figure"] !== nothing
                savefig(p, args["figure"])
                println("Plots saved to: $(args["figure"])")
            else
                display(p)
            end
        catch e
            println("Warning: Could not generate plots. Make sure Plots.jl is installed.")
            println("Error: $e")
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
