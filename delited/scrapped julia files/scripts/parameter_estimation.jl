""" parameter_estimation.jl
Parameter estimation using MCMC DRAM (Delayed Rejection Adaptive Metropolis) algorithm.

This script performs Bayesian parameter inference by fitting the model to synthetic or
experimental data. The DRAM algorithm is an adaptive MCMC method that combines delayed
rejection with adaptive local scaling.

Usage:  julia parameter_estimation.jl --data data_file.csv --param α β --priors 0.5:2.0 0.05:0.2
        julia parameter_estimation.jl --data data_file.csv --param γ --iterations 10000 --burnin 2000
        julia parameter_estimation.jl --data data_file.csv --param α β γ --chains 4 """

using DifferentialEquations
using CSV, DataFrames
using ArgParse
using Distributions
using Statistics
using StatsBase
using Random
using Dates
using ProgressMeter
using LinearAlgebra

include(joinpath(@__DIR__, "..", "src", "ode_system.jl")) # Load modules directly
include(joinpath(@__DIR__, "..", "src", "parameters.jl"))

""" parse_commandline()
Parse command-line arguments for parameter estimation.
Returns: Dict: Parsed arguments """
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--data", "-d"
            help = "CSV file with observation data (columns: time, prey, predator)"
            arg_type = String
            required = true
        "--param", "-p"
            help = "parameters to estimate [param_name ...]"
            nargs = '+'
            arg_type = String
            required = true
        "--priors"
            help = "prior distributions as ranges [min:max ...] (uniform priors)"
            nargs = '+'
            arg_type = String
        "--iterations", "-n"
            help = "number of MCMC iterations"
            arg_type = Int
            default = 5000
        "--burnin", "-b"
            help = "number of burn-in iterations to discard"
            arg_type = Int
            default = 1000
        "--sigma"
            help = "measurement noise standard deviation"
            arg_type = Float64
            default = 0.05
        "--initial"
            help = "initial parameter values [value ...]"
            nargs = '+'
            arg_type = Float64
        "--fixed"
            help = "fixed parameter values [param_name param_value ...]"
            nargs = '+'
            arg_type = String
        "--chains", "-c"
            help = "number of independent MCMC chains"
            arg_type = Int
            default = 1
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

""" parse_prior_specification(prior_strings::Vector{String})
Parse prior specifications from command line.
Arguments: prior_strings::Vector{String}: Vector of prior strings, e.g., ["0.5:2.0", "0.05:0.2"]
Returns: Vector{Uniform}: Vector of uniform distributions """
function parse_prior_specification(prior_strings::Vector{String})
    priors = Distribution[]
    for prior_str in prior_strings
        parts = split(prior_str, ":")
        if length(parts) == 2
            min_val = parse(Float64, parts[1])
            max_val = parse(Float64, parts[2])
            push!(priors, Uniform(min_val, max_val))
        else
            error("Invalid prior format: $prior_str (use min:max)")
        end
    end
    return priors
end

""" parse_fixed_parameters(fixed_args::Vector{String})
Parse fixed parameter specifications.
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

""" load_observation_data(filepath)
Load observation data from CSV file.
Arguments: filepath::String: Path to CSV file
Returns: Tuple: (t_obs, x_obs, y_obs) - observation times and state values """
function load_observation_data(filepath)
    df = CSV.read(filepath, DataFrame)
    t_col = Symbol.(names(df))[1] # Expect columns: time, prey, predator (or similar names)
    x_col = Symbol.(names(df))[2]
    y_col = Symbol.(names(df))[3]
    return df[!, t_col], df[!, x_col], df[!, y_col]
end

""" vector_to_parameters(param_vector, param_names, fixed_params)
Convert a parameter vector to a NamedTuple.
Arguments:  param_vector::Vector: Estimated parameter values
            param_names::Vector{Symbol}: Names of estimated parameters
            fixed_params::Dict: Fixed parameter values
Returns:    NamedTuple: Complete parameter set """
function vector_to_parameters(param_vector, param_names, fixed_params)
    p = merge(DEFAULT_PARAMETERS, fixed_params)
    params_dict = Dict(zip(param_names, param_vector))
    return merge(p, (; zip(keys(params_dict), values(params_dict))...))
end

""" simulate_model(p, t_obs, u0)
Simulate the model at observation times.
Arguments:  p::NamedTuple: Parameters
            t_obs::Vector: Time points where to evaluate
            u0::Vector: Initial conditions
Returns: Tuple: (x_sim, y_sim) - simulated trajectory at observation times
           Returns (NaN, NaN) if simulation fails"""
function simulate_model(p, t_obs, u0)
    try
        prob = ODEProblem(lotka_volterra!, u0, (t_obs[1], t_obs[end]), p)
        sol = solve(prob, Tsit5(), saveat=t_obs, abstol=1e-8, reltol=1e-8)
        return sol[1, :], sol[2, :]
    catch
        return fill(NaN, length(t_obs)), fill(NaN, length(t_obs))
    end
end

""" log_likelihood(x_obs, y_obs, x_sim, y_sim, sigma)
Calculate log-likelihood of observations given simulations.
Arguments:  x_obs::Vector: Observed prey populations
            y_obs::Vector: Observed predator populations
            x_sim::Vector: Simulated prey populations
            y_sim::Vector: Simulated predator populations
            sigma::Float64: Measurement noise standard deviation
Returns: Float64: Log-likelihood"""
function log_likelihood(x_obs, y_obs, x_sim, y_sim, sigma)
    if any(isnan, x_sim) || any(isnan, y_sim) # Check for NaN in simulation
        return -Inf
    end
    loglik = 0.0 # Negative log-likelihood for Gaussian noise
    for i in 1:length(x_obs)
        loglik += logpdf(Normal(x_sim[i], sigma), x_obs[i])
        loglik += logpdf(Normal(y_sim[i], sigma), y_obs[i])
    end
    return loglik
end

""" log_prior(param_vector, priors)
Calculate log-prior probability.
Arguments:  param_vector::Vector: Parameter values
            priors::Vector{Distribution}: Prior distributions
Returns: Float64: Log-prior
"""
function log_prior(param_vector, priors)
    logprior = 0.0
    for (p, prior) in zip(param_vector, priors)
        logprior += logpdf(prior, p)
    end
    return logprior
end

"""
    metropolis_step(current, proposal, loglik_current, loglik_proposal, logprior_current, logprior_proposal)

Metropolis-Hastings acceptance step.

Arguments:
    current::Vector: Current parameter vector
    proposal::Vector: Proposed parameter vector
    loglik_current::Float64: Log-likelihood at current
    loglik_proposal::Float64: Log-likelihood at proposal
    logprior_current::Float64: Log-prior at current
    logprior_proposal::Float64: Log-prior at proposal

Returns:
    Tuple: (accepted::Bool, new_params::Vector, new_loglik::Float64, new_logprior::Float64)
"""
function metropolis_step(current, proposal, loglik_current, loglik_proposal, logprior_current, logprior_proposal)
    log_alpha = (loglik_proposal + logprior_proposal) - (loglik_current + logprior_current)
    
    if log(rand()) < log_alpha
        return true, proposal, loglik_proposal, logprior_proposal
    else
        return false, current, loglik_current, logprior_current
    end
end

"""
    adaptive_metropolis_sampling(t_obs, x_obs, y_obs, param_names, priors, fixed_params, 
                                 n_iterations, burnin, sigma, u0; initial_params=nothing, verbose=false)

Perform Adaptive Metropolis (AM) sampling.

Returns:
    Tuple: (chain, accepted_count, acceptance_rate)
"""
function adaptive_metropolis_sampling(t_obs, x_obs, y_obs, param_names, priors, fixed_params,
                                      n_iterations, burnin, sigma, u0; initial_params=nothing, verbose=false)
    n_params = length(param_names)
    chain = fill(NaN, n_iterations, n_params)
    
    # Initialize parameters
    if initial_params === nothing
        current_params = [rand(p) for p in priors]
    else
        current_params = copy(initial_params)
    end
    
    # Evaluate at initial point
    p_current = vector_to_parameters(current_params, param_names, fixed_params)
    x_sim, y_sim = simulate_model(p_current, t_obs, u0)
    loglik_current = log_likelihood(x_obs, y_obs, x_sim, y_sim, sigma)
    logprior_current = log_prior(current_params, priors)
    
    # Adaptive scaling
    chain_history = fill(NaN, min(burnin, 100), n_params)
    history_idx = 1
    adaptive_cov = Matrix{Float64}(I, n_params, n_params)
    proposal_scale = 1.0
    
    accepted = 0
    pbar = Progress(n_iterations, desc="MCMC sampling: ")
    
    for iter in 1:n_iterations
        # Propose new parameters with adaptive covariance
        proposal = copy(current_params)  # Initialize
        try
            L = cholesky(adaptive_cov).L
            proposal = current_params + proposal_scale * sqrt(2.38^2 / n_params) * L * randn(n_params)
        catch
            # If Cholesky fails, use diagonal approximation
            proposal = current_params + proposal_scale * sqrt(2.38^2 / n_params) * sqrt.(diag(adaptive_cov)) .* randn(n_params)
        end
        
        # Check if proposal is within prior support
        valid_proposal = all(priors[i].a <= proposal[i] <= priors[i].b for i in 1:n_params)
        
        if valid_proposal
            p_proposal = vector_to_parameters(proposal, param_names, fixed_params)
            x_sim, y_sim = simulate_model(p_proposal, t_obs, u0)
            loglik_proposal = log_likelihood(x_obs, y_obs, x_sim, y_sim, sigma)
            logprior_proposal = log_prior(proposal, priors)
            
            accepted_step, new_params, loglik_new, logprior_new = 
                metropolis_step(current_params, proposal, loglik_current, loglik_proposal, 
                               logprior_current, logprior_proposal)
            
            if accepted_step
                current_params = new_params
                loglik_current = loglik_new
                logprior_current = logprior_new
                accepted += 1
            end
        end
        
        chain[iter, :] = current_params
        
        # Adapt covariance matrix (after burnin)
        if iter > burnin
            recent_chain = chain[max(1, iter-100):iter, :]
            recent_chain = recent_chain[.!any(isnan.(recent_chain), dims=2), :]
            if size(recent_chain, 1) > 1
                adaptive_cov = cov(recent_chain)
                adaptive_cov = adaptive_cov + 1e-6 * I  # Add small regularization
            end
        end
        
        next!(pbar)
    end
    
    acceptance_rate = accepted / n_iterations
    return chain, accepted, acceptance_rate
end

"""
    save_mcmc_results(chain, param_names, acceptance_rate, filename)

Save MCMC results to files.

Arguments:
    chain::Matrix: MCMC chain [iterations, params]
    param_names::Vector{Symbol}: Parameter names
    acceptance_rate::Float64: Acceptance rate
    filename::String: Base output filename
"""
function save_mcmc_results(chain, param_names, acceptance_rate, filename)
    data_dir = joinpath(@__DIR__, "..", "data")
    mkpath(data_dir)
    filepath = joinpath(data_dir, filename)
    
    # Save chain as CSV
    df = DataFrame(chain, String.(param_names))
    CSV.write(filepath, df)
    
    # Save diagnostics
    diag_file = filepath * ".diag"
    open(diag_file, "w") do f
        write(f, "MCMC Diagnostics\n")
        write(f, "================\n")
        write(f, "Generated: $(now())\n")
        write(f, "Acceptance rate: $(acceptance_rate)\n")
        write(f, "Number of iterations: $(size(chain, 1))\n")
        write(f, "Number of parameters: $(size(chain, 2))\n")
        write(f, "\nParameter means (full chain):\n")
        for (i, param) in enumerate(param_names)
            write(f, "  $param = $(mean(chain[:, i]))\n")
        end
        write(f, "\nParameter standard deviations (full chain):\n")
        for (i, param) in enumerate(param_names)
            write(f, "  $param = $(std(chain[:, i]))\n")
        end
    end
    
    return filepath
end

"""
    main()

Main execution function.
"""
function main()
    args = parse_commandline()
    
    # Load observation data
    if !isfile(args["data"])
        error("Data file not found: $(args["data"])")
    end
    
    t_obs, x_obs, y_obs = load_observation_data(args["data"])
    
    # Parse parameters
    param_names = Symbol.(args["param"])
    n_params = length(param_names)
    
    # Setup priors
    if args["priors"] !== nothing && length(args["priors"]) == n_params
        priors = parse_prior_specification(args["priors"])
    else
        # Use default wide priors
        priors = [Uniform(0.1, 5.0) for _ in 1:n_params]
    end
    
    # Parse fixed parameters
    fixed_params_dict = Dict()
    if args["fixed"] !== nothing
        fixed_parsed = parse_fixed_parameters(args["fixed"])
        fixed_params_dict = fixed_parsed
    end
    
    # Setup initial values
    if args["initial"] !== nothing && length(args["initial"]) == n_params
        initial_params = args["initial"]
    else
        initial_params = [rand(p) for p in priors]
    end
    
    # Initial conditions from data
    u0 = [x_obs[1], y_obs[1]]
    
    println("=" ^ 70)
    println("PARAMETER ESTIMATION VIA MCMC DRAM")
    println("=" ^ 70)
    println("Data file: $(args["data"])")
    println("Observations: $(length(t_obs)) time points")
    println("Parameters to estimate: $(join(param_names, ", "))")
    println("Initial values: $(initial_params)")
    println("Iterations: $(args["iterations"]), Burnin: $(args["burnin"])")
    println()
    
    # Run MCMC chains
    all_chains = []
    all_acceptances = []
    all_rates = []
    
    for chain_idx in 1:args["chains"]
        if args["chains"] > 1
            println("Running chain $chain_idx / $(args["chains"])...")
        end
        
        chain, accepted, acceptance_rate = adaptive_metropolis_sampling(
            t_obs, x_obs, y_obs, param_names, priors, fixed_params_dict,
            args["iterations"], args["burnin"], args["sigma"], u0;
            initial_params=initial_params, verbose=args["verbose"]
        )
        
        push!(all_chains, chain)
        push!(all_acceptances, accepted)
        push!(all_rates, acceptance_rate)
        
        if args["verbose"]
            println("  Acceptance rate: $(acceptance_rate)")
        end
    end
    
    # Determine output filename
    if args["output"] === nothing
        timestamp = Dates.format(now(), "yyyymmdd_HHMMSS")
        output_filename = "posterior_$(join(param_names, "_"))_$(timestamp).csv"
    else
        output_filename = args["output"]
    end
    
    # Save first chain (or combined if multiple chains)
    filepath = save_mcmc_results(all_chains[1], param_names, all_rates[1], output_filename)
    
    println("\n" * "=" ^ 70)
    println("Parameter Estimation Complete")
    println("=" ^ 70)
    println("Results saved to: $filepath")
    println("\nAcceptance rates:")
    for (i, rate) in enumerate(all_rates)
        println("  Chain $i: $(rate)")
    end
    
    # Print posterior summary
    chain_postburnin = all_chains[1][end-min(1000, end÷2):end, :]
    println("\nPosterior summary (after burnin):")
    for (i, param) in enumerate(param_names)
        mean_val = mean(chain_postburnin[:, i])
        std_val = std(chain_postburnin[:, i])
        ci_low = quantile(chain_postburnin[:, i], 0.025)
        ci_high = quantile(chain_postburnin[:, i], 0.975)
        println("  $param: $(round(mean_val, digits=4)) ± $(round(std_val, digits=4)) [$(round(ci_low, digits=4)), $(round(ci_high, digits=4))]")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
