""" Simulation Script for Lotka-Volterra System
Usage:
    julia simulate.jl --help                  # Show available options
    julia simulate.jl --time 100              # Simulate for 100 time units
    julia simulate.jl --settle --x0 1 2 --p default  # Run until equilibrium
    julia simulate.jl --settle --x0 1 2 --α 1.2 --β 0.1  # Custom parameters
Examples of output files are saved in ../data/ """

using DifferentialEquations
using CSV, DataFrames
using ArgParse
using Statistics
using Dates

include(joinpath(@__DIR__, "..", "src", "TIME_ode_model.jl")) # Load modules directly
include(joinpath(@__DIR__, "..", "src", "TIME_params.jl"))

""" parse_commandline()
Parse command-line arguments for the simulation script.
Returns: Dict: Parsed arguments
"""
function parse_commandline()
    s = ArgParseSettings() 
    @add_arg_table! s begin
        "--time", "-t"
            help = "simulation time (run for fixed duration)"
            arg_type = Float64
            default = nothing
        "--settle", "-s"
            help = "run until equilibrium/limit cycle is detected"
            action = :store_true
        "--x0"
            help = "initial conditions [x[i] for i in list of state variables]"
            nargs = 2
            arg_type = Float64
            default = [1.0, 0, 0, 0, 85000, 15000, 71000, 12000, 56000, 8000, 0.0085, 0.12, 0.0094]
        "--params", "-p"
            help = "predefined parameter set name (default, fast_predator, strong_predation, weak_predation)"
            arg_type = String
            default = "default" 
        "--output", "-o"
            help = "output filename (in data/ directory)"
            arg_type = String
            default = nothing
        "--verbose", "-v"
            help = "print detailed output"
            action = :store_true
    end
    return parse_args(s)
end

""" setup_parameters(args::Dict)
Setup parameter set from command-line arguments.
Arguments: args::Dict: Parsed command-line arguments
Returns: NamedTuple: Final parameter set
"""
function setup_parameters(args::Dict)
    params = get_parameter_set(args["params"])
    override_params = Dict() # Override with command-line specified parameters
    for param in (:beta_M2, :beta_Tc, :beta_Th1CK2, :beta_Th1CK3, :beta_Th2, :beta_Treg, :gamma_C, :gamma_CR, :gamma_M1, :gamma_M2, :gamma_S,  :gamma_Tc, :gamma_Th1, :gamma_Th2, :gamma_Treg, :delta_C, :delta_Ck1, :delta_Ck2, :delta_Ck3, :delta_CR, :delta_M1, :delta_M2, :delta_S, :delta_Tc, :delta_Th1, :delta_Th2, :delta_Treg, :lambda_M1, :lambda_M2, :lambda_Tc1, :lambda_Tc2, :lambda_Tc3, :lambda_Tc4, :lambda_Th1, :lambda_Th2, :lambda_Treg2, :mu_C1, :mu_C2, :mu_S, :mu_SR, :mu_TcS, :mu_TcTreg, :mu_Th1Ck1, :mu_Th1Ck3, :mu_TregCk1, :C_max, :CR_max, :k1, :k11, :k2, :k3, :k4, :k5, :k6, :k8, :k9, :ktc1, :ktc2, :ktc3, :ktc4, :m_C, :m_S, :p_1, :p_2, :r_1, :r_2, :tck, :mu_M1Ck2, :mu_M2Ck1, :k_7, :k_10)
        if haskey(args, string(param)) && args[string(param)] !== nothing
            override_params[param] = args[string(param)]
        end
    end
    if !isempty(override_params)
        params = merge(params, (; zip(keys(override_params), values(override_params))...))
    end
    return params
end

""" convergence_callback(transient_threshold=10.0, check_period=10.0)
Create a callback function to detect convergence to equilibrium/limit cycle.
Arguments: transient_threshold::Float64: Time to wait before checking convergence
            check_period::Float64: Time window to check for convergence
Returns: DiscreteCallback: Callback for ODE solver
"""
function convergence_callback(transient_threshold=10.0, check_period=10.0)
    history = Dict("times" => Float64[], "states" => SVector[]) # Store history of recent states
    
    function condition(u, t, integrator)
        if t > transient_threshold
            push!(history["times"], t)
            push!(history["states"], u)
            if length(history["times"]) > 2 # Keep only last check_period worth of states
                if history["times"][end] - history["times"][1] > check_period
                    popfirst!(history["times"])
                    popfirst!(history["states"])
                end
            end
            if length(history["states"]) >= 2 # Check for convergence
                state_diffs = [norm(history["states"][i] - history["states"][1])  # Check if all recent states are similar (within tolerance)
                              for i in 1:length(history["states"])]
                if all(diff < 1e-3 for diff in state_diffs)
                    return true
                end
            end
        end
        return false
    end
    function affect!(integrator)
        terminate!(integrator)
    end
    return DiscreteCallback(condition, affect!, save_positions=(true, true))
end

""" run_simulation_fixed_time(u0, p, tspan; solver_args...)
Run simulation for a fixed time period.
Arguments:  u0::Vector: Initial conditions
            p::NamedTuple: Parameters
            tspan::Tuple: Time span (t_start, t_end)
            solver_args: Additional arguments for ODE solver
Returns: ODESolution: Solution object"""
function run_simulation_fixed_time(u0, p, tspan; solver_args...)
    prob = ODEProblem(TIME_model_1!, u0, tspan, p)
    sol = solve(prob, Tsit5(); saveat=0.1, solver_args...)
    return sol
end

""" run_simulation_until_settle(u0, p, max_time=1000.0; transient=10.0)
Run simulation until system reaches equilibrium or limit cycle.
Arguments:  u0::Vector: Initial conditions
            p::NamedTuple: Parameters
            max_time::Float64: Maximum simulation time
            transient::Float64: Time to wait before checking convergence
Returns: ODESolution: Solution object
"""
function run_simulation_until_settle(u0, p, max_time=1000.0; transient=10.0)
    cb = convergence_callback(transient, 5.0)
    prob = ODEProblem(TIME_model_1!, u0, (0.0, max_time), p)
    sol = solve(prob, Tsit5(); callback=cb, saveat=0.1)
    return sol
end

""" save_solution(sol, p, filename::String)
Save solution to CSV file with metadata.
Arguments:  sol::ODESolution: Solution object
            p::NamedTuple: Parameters used
            filename::String: Output filename (without directory) """
function save_solution(sol, p, filename::String)
    data_dir = joinpath(@__DIR__, "..", "data")
    mkpath(data_dir)
    filepath = joinpath(data_dir, filename)
    df = DataFrame( # Create DataFrame
        time = vec(sol.t), prey = vec(sol[1, :]), predator = vec(sol[2, :]))
    CSV.write(filepath, df) # Save to CSV
    metadata_file = filepath * ".meta" # Also save metadata
    open(metadata_file, "w") do f
        write(f, "Simulation metadata\n")
        write(f, "==================\n")
        write(f, "Generated: $(now())\n")
        write(f, "Final time: $(sol.t[end])\n")
        write(f, "Number of time points: $(length(sol.t))\n")
        write(f, "\nParameters:\n")
        for (key, val) in pairs(p)
            write(f, "  $key = $val\n")
        end
        #write(f, "\nInitial conditions:\n")
        #write(f, "  u (state vector) = $(sol[:, 1])\n")
        #write(f, "\nFinal state:\n")
        #write(f, "  u_final = $(sol[:, end])\n") # check if the syntax is correct for [:,]
    end
    return filepath
end

""" main()
Main execution function. """
function main()
    args = parse_commandline()    
    params = setup_parameters(args) # Setup parameters
    u0 = args["x0"]
    if args["verbose"]
        println("=" ^ 60)
        println("Tumour-Immune Microenvironment Dynamics Simulation")
        println("=" ^ 60)
        println("Initial conditions: state_vector = $(u0)")
        println("Parameters:")
        print_parameter_set(params)
        println()
    end
    if args["time"] !== nothing # Run simulation
        if args["verbose"]
            println("Running simulation for t ∈ [0, $(args["time"])]...")
        end
        sol = run_simulation_fixed_time(u0, params, (0.0, args["time"]))
    elseif args["settle"]
        if args["verbose"]
            println("Running simulation until equilibrium/limit cycle...")
        end
        sol = run_simulation_until_settle(u0, params)
    else
        if args["verbose"]
            println("No simulation mode specified. Running for t ∈ [0, 100]...")
        end
        sol = run_simulation_fixed_time(u0, params, (0.0, 100.0))
    end
    if args["output"] === nothing # Determine output filename
        timestamp = Dates.format(now(), "yyyymmdd_HHMMSS")
        mode = args["settle"] ? "settle" : "fixed"
        output_filename = "trajectory_$(mode)_$(timestamp).csv"
    else
        output_filename = args["output"]
    end
    filepath = save_solution(sol, params, output_filename) # Save solution
    if args["verbose"]
        println("Simulation completed!")
        println("  Final time: $(sol.t[end])")
        #println("  Final state: x = $(sol[:, end]), y = $(sol[2, end])")
        println("  Data saved to: $filepath")
        println()
        println("Try: julia analyze_trajectory.jl --file $output_filename")
    else
        println("Data saved to: $filepath")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
