
using GlobalSensitivity, QuasiMonteCarlo
using DifferentialEquations
include("model.jl")
include("parameters.jl")

using .CancerModel
using .ModelParameters

function run_efast_analysis()
    p_dict = get_default_parameters()
    u0 = get_initial_conditions()
    
    # Define search ranges (±20% of baseline)
    lb = [p_dict[n] * 0.8 for n in param_names]
    ub = [p_dict[n] * 1.2 for n in param_names]
    
    # Objective function: Total tumor mass at day 800
    function f(p_matrix)
        results = zeros(size(p_matrix, 2))
        for i in 1:size(p_matrix, 2)
            # Temp param dict
            p_local = copy(p_dict)
            for (j, name) in enumerate(param_names)
                p_local[name] = p_matrix[j, i]
            end
            
            prob = ODEProblem(tumor_ode!, u0, (0.0, 800.0), p_local)
            sol = solve(prob, saveat=[800.0])
            
            if sol.retcode == :Success
                results[i] = sol.u[end][1] + sol.u[end][2] + sol.u[end][3] + sol.u[end][4]
            else
                results[i] = 1e20 # Penalty for divergence
            end
        end
        return results
    end
    
    # Run eFAST: 100 samples per curve, 5 resamples
    res = gsa(f, eFAST(num_resamples=5), [[lb[i], ub[i]] for i in 1:length(lb)], n=100)
    
    return res
end
