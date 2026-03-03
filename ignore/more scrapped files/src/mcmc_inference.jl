
using AdaptiveMCMC
using Statistics
using CSV
using DataFrames
include("model.jl")
include("parameters.jl")

using .CancerModel
using .ModelParameters

function run_inference(data_csv_path, iterations=500000)
    # Load experimental data
    df = CSV.read(data_csv_path, DataFrame)
    
    # Log-likelihood function
    function log_likelihood(params_vec)
        # Map vec to Dict
        p = get_default_parameters()
        for (i, val) in enumerate(params_vec)
            p[param_names[i]] = val
        end
        
        u0 = get_initial_conditions()
        prob = ODEProblem(tumor_ode!, u0, (0.0, 7.0), p)
        sol = solve(prob, saveat=1.0)
        
        if sol.retcode != :Success
            return -Inf
        end
        
        # Compare C population (u[3]) to Gastric data
        ll = 0.0
        for i in 1:nrow(df)
            t_idx = Int(df.time[i]) + 1
            if t_idx <= length(sol.t)
                ll += -0.5 * (sol.u[t_idx][3] - df.cell_count[i])^2 / 1e10
            end
        end
        return ll
    end

    # Prior: Normal centered at default params
    p_baseline = [get_default_parameters()[n] for n in param_names[1:5]] # Estimate first 5 as example
    
    # Run Adaptive MCMC
    out = adaptive_rwm(p_baseline, log_likelihood, iterations; algorithm=:am)
    
    return out
end
