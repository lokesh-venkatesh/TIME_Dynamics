
module ModelParameters

export get_default_parameters, get_initial_conditions, param_names

# Parameter names corresponding to the 71 variables in the paper
const param_names = [
    "gamma_s", "m_s", "p1", "p2", "delta_s", "mu_s", "k1", "tck", "ktc1", # S dynamics
    "mu_sr", "k2", "ktc2", # Sr dynamics
    "gamma_c", "m_c", "Ktumor", "r1", "delta_c", "mu_c1", "k3", "mu_c2", "k4", "ktc3", # C dynamics
    "r2", "mu_cr", "mu_c1_r", "k5", "mu_c2_r", "k6", "ktc4", # Cr dynamics
    "gamma_m1", "lambda_m1", "delta_m1", "mu_m1ck2", "k7", # M1 dynamics
    "gamma_m2", "lambda_m2", "delta_m2", "mu_m2ck1", "k10", # M2 dynamics
    "gamma_th1", "lambda_th1", "delta_th1", "mu_th1ck1", "k8", "mu_th1ck3", "k9", # Th1
    "gamma_th2", "lambda_th2", "delta_th2", # Th2
    "gamma_tc", "lambda_tc1", "lambda_tc4", "mu_tcs", "lambda_tc2", "delta_tc", "mu_tctreg", "lambda_tc3", # Tc
    "gamma_treg", "lambda_treg2", "delta_treg", "mu_tregck1", "k11", # Treg
    "beta_m2", "delta_ck1", "beta_treg", "beta_th2", # IL10
    "beta_th1ck2", "beta_tc", "delta_ck2", # IFNg
    "beta_th1ck3", "delta_ck3" # IL2
]

function get_default_parameters()
    # Baseline parameters from the paper (representative values)
    p = Dict{String, Float64}()
    
    # Growth & Mutation
    p["gamma_s"] = 0.15
    p["m_s"] = 2.0e-7
    p["p1"] = 0.3
    p["p2"] = 0.1
    p["delta_s"] = 0.02
    
    # Tumor Capacity
    p["Ktumor"] = 1.0e10
    p["r1"] = 10.0
    p["r2"] = 10.0
    
    # Killing / Interactions
    p["mu_s"] = 0.5
    p["k1"] = 1.0e6
    p["tck"] = 0.1
    p["ktc1"] = 1.0e5
    p["ktc2"] = 1.0e5
    p["ktc3"] = 1.0e5
    p["ktc4"] = 1.0e5
    
    # Immune components (fill remaining 71 based on steady state calibration)
    # Note: In a real scenario, these would be loaded from a CSV or S1 Text
    for name in param_names
        if !haskey(p, name)
            p[name] = 0.1 # Placeholder for calibrated values
        end
    end
    
    return p
end

function get_initial_conditions()
    # 13 dimensions: S, Sr, C, Cr, M1, M2, Th1, Th2, Tc, Treg, IL10, IFNg, IL2
    u0 = zeros(13)
    u0[1] = 1.0    # S (One Cancer Stem Cell)
    # Others initialized to 0 or healthy baseline
    u0[5:13] .= [100.0, 100.0, 50.0, 50.0, 100.0, 50.0, 0.01, 0.01, 0.01] 
    return u0
end

end
