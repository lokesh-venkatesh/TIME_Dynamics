""" parameters.jl
This script is supposed to define 
parameter sets and related utility functions"""

""" DEFAULT_PARAMETERS
Default parameter set that may be imported


NOTE: Estimated  parameter  values  have  been  determined  by  the  MCMC  techniques  using  the  time 
course  experiment  cytometric  data  for  cancer  cell  proliferation  for  Gastric  cancer  cell  line 
(SGC7901) (Refer Methods of Main Article) [15]. 
NOTE: Expected  parameter  values  are  estimated  by  varying  the  parameters  within  the  biologically 
feasible ranges found in various Literatures so as to determine its expected value to calibrate the 
model with experimental observations.
"""
const DEFAULT_PARAMETERS = (
    beta_M2 = 1e-15, # Expected
    beta_Tc = 1e-8, # Expected
    beta_Th1CK2 = 1e-7, # Expected
    beta_Th1CK3 = 1e-8, # Expected
    beta_Th2 = 1e-9, # Expected
    beta_Treg = 1e-10, # Expected
    gamma_C = 0.1282, # ESTIMATED
    gamma_CR = 0.1282, # Expected
    gamma_M1 = 0.7, # Expected
    gamma_M2 = 0.01, # Expected
    gamma_S = 0.15, # Expected
    gamma_Tc = 1.0, 
    gamma_Th1 = 2.0, 
    gamma_Th2 = 2.0, 
    gamma_Treg = 0.3, 
    delta_C = 0.8055, # ESTIMATED
    delta_Ck1 = 19.757, # ESTIMATED
    delta_Ck2 = 6.1212, # ESTIMATED
    delta_Ck3 = 8.664339, # CALCULATED
    delta_CR = 5.37e-5, # ESTIMATED
    delta_M1 = 1.02, 
    delta_M2 = 0.05, 
    delta_S = 2e-7, # Expected
    delta_Tc = 5.2939, # ESTIMATED
    delta_Th1 = 2.0, 
    delta_Th2 = 2.0, # Expected
    delta_Treg = 1.0, 
    lambda_M1 = 1e8, # Expected 
    lambda_M2 = 1e6, # Expected 
    lambda_Tc1 = 1e5, # Expected
    lambda_Tc2 = 5e5, # Expected
    lambda_Tc3 = 5e10, # Expected
    lambda_Tc4 = 1e5, # Expected
    lambda_Th1 = 1e5, # Expected
    lambda_Th2 = 1e5, # Expected
    lambda_Treg2 = 1e7, # Expected
    mu_C1 = 0.75, 
    mu_C2 = 0.9, 
    mu_S = 0.17, 
    mu_SR = 0.18, # Expected
    mu_TcS = 1e-10, 
    mu_TcTreg = 1.5e-5, 
    mu_Th1Ck1 = 1e-9, 
    mu_Th1Ck3 = 0.1245, 
    mu_TregCk1 = 1e-7, # Expected
    C_max = 1e10, 
    CR_max = 1e10, # Expected
    k1 = 10.0, # Expected
    k11 = 0.001, # Expected
    k2 = 10.0, # Expected
    k3 = 2.0531, # ESTIMATED
    k4 = 3.02, # ESTIMATED
    k5 = 6.7979, # ESTIMATED
    k6 = 6.9937, # ESTIMATED
    k8 = 0.01, # Expected
    k9 = 0.001, # Expected
    ktc1 = 1e9, # Expected
    ktc2 = 1e8, # Expected
    ktc3 = 1e9, # Expected
    ktc4 = 1e9, # Expected
    m_C = 0.01, # Expected
    m_S = 4e-7, 
    p_1 = 0.2, 
    p_2 = 0.05, 
    r_1 = 0.0001, # Expected
    r_2 = 1e-5, # Expected
    tck = 0.1, # Expected
    mu_M1Ck2 = 0.01, # Expected
    mu_M2Ck1 = 0.02, # Expected
    k_7 = 0.2, # Expected
    k_10 = 0.0 # Expected
)


""" get_parameter_set()
Retrieve a pre-defined parameter set by name."""
function get_parameter_set(name::String)
    sets = Dict(
        "default" => DEFAULT_PARAMETERS,
    )    
    if !haskey(sets, name)
        error("Unknown parameter set: $name. Available: $(join(keys(sets), ", "))")
    end
    return sets[name]
end


""" create_parameter_set()
Create a custom parameter set (with the default being DEFAULT_PARAMETERS)"""
function create_parameter_set(; 
    beta_M2 = 1e-15, # Expected
    beta_Tc = 1e-8, # Expected
    beta_Th1CK2 = 1e-7, # Expected
    beta_Th1CK3 = 1e-8, # Expected
    beta_Th2 = 1e-9, # Expected
    beta_Treg = 1e-10, # Expected
    gamma_C = 0.1282, # ESTIMATED
    gamma_CR = 0.1282, # Expected
    gamma_M1 = 0.7, # Expected
    gamma_M2 = 0.01, # Expected
    gamma_S = 0.15, # Expected
    gamma_Tc = 1.0, 
    gamma_Th1 = 2.0, 
    gamma_Th2 = 2.0, 
    gamma_Treg = 0.3, 
    delta_C = 0.8055, # ESTIMATED
    delta_Ck1 = 19.757, # ESTIMATED
    delta_Ck2 = 6.1212, # ESTIMATED
    delta_Ck3 = 8.664339, # CALCULATED
    delta_CR = 5.37e-5, # ESTIMATED
    delta_M1 = 1.02, 
    delta_M2 = 0.05, 
    delta_S = 2e-7, # Expected
    delta_Tc = 5.2939, # ESTIMATED
    delta_Th1 = 2.0, 
    delta_Th2 = 2.0, # Expected
    delta_Treg = 1.0, 
    lambda_M1 = 1e8, # Expected 
    lambda_M2 = 1e6, # Expected 
    lambda_Tc1 = 1e5, # Expected
    lambda_Tc2 = 5e5, # Expected
    lambda_Tc3 = 5e10, # Expected
    lambda_Tc4 = 1e5, # Expected
    lambda_Th1 = 1e5, # Expected
    lambda_Th2 = 1e5, # Expected
    lambda_Treg2 = 1e7, # Expected
    mu_C1 = 0.75, 
    mu_C2 = 0.9, 
    mu_S = 0.17, 
    mu_SR = 0.18, # Expected
    mu_TcS = 1e-10, 
    mu_TcTreg = 1.5e-5, 
    mu_Th1Ck1 = 1e-9, 
    mu_Th1Ck3 = 0.1245, 
    my_TregCk1 = 1e-7, # Expected
    C_max = 1e10, 
    CR_max = 1e10, # Expected
    k1 = 10.0, # Expected
    k11 = 0.001, # Expected
    k2 = 10.0, # Expected
    k3 = 2.0531, # ESTIMATED
    k4 = 3.02, # ESTIMATED
    k5 = 6.7979, # ESTIMATED
    k6 = 6.9937, # ESTIMATED
    k8 = 0.01, # Expected
    k9 = 0.001, # Expected
    ktc1 = 1e9, # Expected
    ktc2 = 1e8, # Expected
    ktc3 = 1e9, # Expected
    ktc4 = 1e9, # Expected
    m_C = 0.01, # Expected
    m_S = 4e-7, 
    p_1 = 0.2, 
    p_2 = 0.05, 
    r_1 = 0.0001, # Expected
    r_2 = 1e-5, # Expected
    tck = 0.1, # Expected
    mu_M1Ck2 = 0.01, # Expected
    mu_M2Ck1 = 0.02, # Expected
    k_7 = 0.2, # Expected
    k_10 = 0.0 # Expected
)
    return (beta_M2 = beta_M2, 
    beta_Tc = beta_Tc,
    beta_Th1CK2 = beta_Th1CK2,
    beta_Th1CK3 = beta_Th1CK3,
    beta_Th2 = beta_Th2 ,
    beta_Treg = beta_Treg, 
    gamma_C = gamma_C, 
    gamma_CR = gamma_CR, 
    gamma_M1 = gamma_M1, 
    gamma_M2 = gamma_M2,
    gamma_S = gamma_S, 
    gamma_Tc = gamma_Tc,
    gamma_Th1 = gamma_Th1,
    gamma_Th2 = gamma_Th2,
    gamma_Treg = gamma_Treg,
    delta_C = delta_C,
    delta_Ck1 = delta_Ck1,
    delta_Ck2 = delta_Ck2,
    delta_Ck3 = delta_Ck3,
    delta_CR = delta_CR,
    delta_M1 = delta_M1,
    delta_M2 = delta_M2,
    delta_S = delta_S,
    delta_Tc = delta_Tc,
    delta_Th1 = delta_Th1,
    delta_Th2 = delta_Th2,
    delta_Treg = delta_Treg,
    lambda_M1 = lambda_M1,
    lambda_M2 = lambda_M2,
    lambda_Tc1 = lambda_Tc1,
    lambda_Tc2 = lambda_Tc2,
    lambda_Tc3 = lambda_Tc3,
    lambda_Tc4 = lambda_Tc4,
    lambda_Th1 = lambda_Th1,
    lambda_Th2 = lambda_Th2,
    lambda_Treg2 = lambda_Treg2,
    mu_C1 = mu_C1, 
    mu_C2 = mu_C2, 
    mu_S = mu_S, 
    mu_SR = mu_SR, # Expected
    mu_TcS = mu_TcS, 
    mu_TcTreg = mu_TcTreg,
    mu_Th1Ck1 = mu_Th1Ck1, 
    mu_Th1Ck3 = mu_Th1Ck1, 
    mu_TregCk1 = mu_TcTregCk1, # Expected
    C_max = C_max, 
    CR_max = CR_max, # Expected
    k1 = k1, # Expected
    k11 = k11, # Expected
    k2 = k2, # Expected
    k3 = k3, # ESTIMATED
    k4 = k4, # ESTIMATED
    k5 = k5, # ESTIMATED
    k6 = k6, # ESTIMATED
    k8 = k8, # Expected
    k9 = k9, # Expected
    ktc1 = ktc1, # Expected
    ktc2 = kct2, # Expected
    ktc3 = kct3, # Expected
    ktc4 = kct4, # Expected
    m_C = m_C, # Expected
    m_S = m_S, 
    p_1 = p_1, 
    p_2 = p_2, 
    r_1 = r_1, # Expected
    r_2 = r_2, # Expected
    tck = tck, # Expected
    mu_M1Ck2 = mu_M1Ck2, # Expected
    mu_M2Ck1 = mu_M2Ck1, # Expected
    k_7 = k_7, # Expected
    k_10 = k_10 # Expected
    )
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