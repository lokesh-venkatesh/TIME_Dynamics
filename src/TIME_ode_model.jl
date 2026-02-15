""" TIME_ode_model.jl (<this> entire script)
Supposed to contain the model equations for the ODE system under analysis."""

""" TIME_model_1!(du, u, p, t)
This is for updating the ODE system in-place
    
du are the time-derivatives modified in place
u is the vector of state variables
p is the tuple of parameters defined beforehand
t is time (implicit in this particular ODE system)"""
function TIME_model_1!(dx, x, params, t)
    # declaring state variables
    S, S_R, C, C_R, M_1, M_2, T_H1, T_H2, T_C, T_reg, IL10, IFN_gamma, IL2 = x
    # declaring parameters
    beta_M2, beta_Tc, beta_Th1CK2, beta_Th1CK3, beta_Th2, beta_Treg, gamma_C, gamma_CR, gamma_M1, gamma_M2, gamma_S, gamma_Tc, gamma_Th1, gamma_Th2, gamma_Treg, delta_C, delta_Ck1, delta_Ck2, delta_Ck3, delta_CR, delta_M1, delta_M2, delta_S, delta_Tc, delta_Th1, delta_Th2, delta_Treg, lambda_M1, lambda_M2, lambda_Tc1, lambda_Tc2, lambda_Tc3, lambda_Tc4, lambda_Th1, lambda_Th2, lambda_Treg2, mu_C1, mu_C2, mu_S, mu_SR, mu_TcS, mu_TcTreg, mu_Th1Ck1, mu_Th1Ck3, mu_TregCk1, C_max, CR_max, k1, k11, k2, k3, k4, k5, k6, k8, k9, ktc1, ktc2, ktc3, ktc4, m_C, m_S, p_1, p_2, r_1, r_2, tck, mu_M1Ck2, mu_M2Ck1, k_7, k_10   =   params.beta_M2, params.beta_Tc, params.beta_Th1CK2, params.beta_Th1CK3, params.beta_Th2, params.beta_Treg, params.gamma_C, params.gamma_CR, params.gamma_M1, params.gamma_M2, params.gamma_S, params.gamma_Tc, params.gamma_Th1, params.gamma_Th2, params.gamma_Treg, params.delta_C, params.delta_Ck1, params.delta_Ck2, params.delta_Ck3, params.delta_CR, params.delta_M1, params.delta_M2, params.delta_S, params.delta_Tc, params.delta_Th1, params.delta_Th2, params.delta_Treg, params.lambda_M1, params.lambda_M2, params.lambda_Tc1, params.lambda_Tc2, params.lambda_Tc3, params.lambda_Tc4, params.p.lambda_Th1, params.lambda_Th2, params.lambda_Treg2, params.mu_C1, params.mu_C2, params.mu_S, params.mu_SR, params.mu_TcS, params.mu_TcTreg, params.mu_Th1Ck1, params.mu_Th1Ck3, params.mu_TregCK1, params.C_max, params.CR_max, params.k1, params.k11, params.k2, params.k3, params.k4, params.k5, params.k6, params.k8, params.k9, params.ktc1, params.ktc2, params.ktc3, params.ktc4, params.m_C, params.m_S, params.p_1, params.p_2, params.r_1, params.r_2, params.tck, params.mu_M1Ck2, params.mu_M2Ck1, params.k_7, params.k_10
    # governing equations for the dynamical system
    # for the four tumour cell populations
    dx[1] = (gamma_S*(1-m_S)*(1-p_1-p_2))*S - (delta_S+(p_2*gamma_S)+gamma_S*m_S*p_1/2)*S - (mu_S*S*IFN_gamma)/(k1 + IFN_gamma) - (tck*S*T_C)/(ktc1+T_C)
    dx[2] = (gamma_S*(1-p_1-p_2) - (delta_S+(p_2*gamma_S)))*S_R + m_S*gamma_S*(1-p_1/2-p_2)*S - (mu_SR*S_R*IFN_gamma)/(k2 + IFN_gamma) - (tck*S_R*T_C)/(ktc2+T_C)
    dx[3] = gamma_C*(1-m_C)*log((C_max)/(C+r_1))*C + gamma_S*(p_1+p_2)*S - delta_C*C - m_C*gamma_C*C + (mu_C1*C*IL10)/(IL10+k3) - (mu_C2*C*IFN_gamma)/(IFN_gamma+k4) - (tck*C*T_C)/(ktc3+T_C)
    dx[4] = gamma_C*C_R*log((CR_max)/(C_R+r_2)) + gamma_S*S_R*(p_1+p_2) + m_C*gamma_C*C - delta_CR*C_R + (mu_C1*C_R*IL10)/(IL10+k5) - (mu_C2*C_R*IFN_gamma)/(IFN_gamma+k6) - (tck*C_R*T_C)/(ktc4+T_C)
    # NOTE: 0.5*K_tumor = C_max = CR_max
    # for the effector populations:
    dx[5] = gamma_M1*M_1*((C+C_R)/(M_1+lambda_M1)) - delta_M1*M_1 + ((mu_M1Ck2*M1*IFN_gamma)/(IFN_gamma+k_7))
    dx[6] = gamma_M2*M_2*((C+C_R)/(M_2+lambda_M2)) - delta_M2*M_2 + ((mu_M2Ck1*M2*IL10)/(IL10+k_10))
    dx[7] = gamma_Th1*((T_H1*M_1)/(lambda_Th1+T_H1)) - delta_Th1*T_H1 - (mu_Th1Ck1*IL10*T_H1)/(IL10+k8) + (mu_Th1Ck3*IL2*T_H1)/(IL2+k9)
    dx[8] = gamma_Th2*((T_H2*M_2)/(lambda_Th2+T_H2)) - delta_Th2*T_H2
    dx[9] = gamma_Tc*T_C*((C+C_R)/(T_C+lambda_Tc1)) + gamma_Tc*((T_C*T_H1)/(T_C+lambda_Tc4)) - mu_TcS*T_C*((S+S_R)/(T_C+lambda_Tc2)) - delta_Tc*T_C - mu_TcTreg*T_C*((T_reg)/(lambda_Tc3+T_reg))
    # for the inhibitor populations:
    dx[10] = gamma_Treg*((T_reg*M_2)/(T_reg+lambda_Treg2)) - delta_Treg*T_reg + mu_TregCK1*((IL10*T_reg)/(T_reg+k11))
    dx[11] = beta_M2*M2 - delta_Ck1*IL10 + beta_Treg*T_reg + beta_Th2*T_H2
    dx[12] = beta_Th1CK2*T_H1 + beta_Tc*T_C - delta_Ck2*IFN_gamma
    dx[13] = beta_Th1CK3*T_H1 - delta_Ck3*IL2
end


""" TIME_model_1(du, u, p, t)
Non-mutating version of the Lotka-Volterra system.
Returns: Vector of time derivatives [dx/dt, dy/dt] """
function TIME_model_1(x, params, t)
    dx = similar(x)
    TIME_model_1!(dx, x, params, t)
    return dx
end


""" state_has_equilibrium(x, params; eps=1e-6)
Checks if a given state vector is an equilibrium
Arguments:  x::Vector: State variables
        params::NamedTuple: Parameters
        tolerance::Float64: Tolerance for checking if derivatives are near zero
Returns:    Bool: true if the state is near an equilibrium """
function state_has_equilibrium(x, params; tolerance=1e-6)
    dx = similar(x)
    TIME_model_1!(dx, x, params, 0.0)
    return all(abs.(dx) .< tolerance)
end


""" get_state_summary(x)
Provide a human-readable summary of the current state.
Arguments: u::Vector: State variables [x, y]
Returns: String: Formatted summary """
function get_state_summary(u)
    return "S (Cancer Stem Cells): $(u[1]), S_R (Resistant Cancer Stem Cells): $(u[2]), C (Cancer Cells): $(u[3]), C_R (Resistant Cancer Cells): $(u[4]), M_1 (Type-I Tumour Associated Macrophages): $(u[5]), M_2 (Type-II Tumour Associated Macrophages): $(u[6]), T_H1 (Type-I Helper T-Cell): $(u[7]), T_H2 (Type-II Helper T-Cell): $(u[8]), T_C (Cytotoxic T-Cell): $(u[9]), T_reg (Regulatory T-Cell): $(u[10]), IL10 (Interleukin-10): $(u[11]), IFN_gamma (Interferon-gamma): $(u[12]), IL2 (Interleukin-2): $(u[13])"
end

export TIME_model_1!, TIME_model_1, state_has_equilibrium, get_state_summary