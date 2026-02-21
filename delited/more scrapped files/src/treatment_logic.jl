
module TreatmentProtocols

using DifferentialEquations

export apply_protocol_1, apply_protocol_2

# Therapy intensity functions
# u1: Radiotherapy probability (Eq 14)
# u2_S, u2_C: Chemotherapy kill probability (Eq 15-16)
# u3_Tc, u3_Th1: Immunotherapy boost (Eq 17-18)

function get_radiotherapy_effect(dose)
    alpha = 0.01; beta = 0.001 # Radiosensitivity params
    return 1.0 - exp(-alpha * dose - beta * dose^2)
end

function get_chemo_effect_S(dose, freq, ks)
    M = 0.5; Mc = 0.1 # Drug efficiency
    return freq * (1.0 - exp(-Mc * dose)) - ks
end

function get_chemo_effect_C(dose, freq)
    Mc = 0.1
    return freq * (1.0 - exp(-Mc * dose))
end

function apply_protocol_1(prob)
    # Protocol 1: DT200 -> Ch(14d, 6 cycles) -> R(28d) -> FT(15d) -> Ch(14d, 6 cycles)
    # This involves discrete changes to the state variables (cell death events)
    # We use a CallbackSet to manage these pulses.
    
    callbacks = []
    
    # Chemotherapy cycles (6 cycles, 14 days each) starting at day 200
    for cycle in 0:5
        start_t = 200.0 + cycle * 14.0
        cb = PresetTimeCallback(start_t, (integrator) -> begin
            u2_S = get_chemo_effect_S(800, 1, 0.01)
            u2_C = get_chemo_effect_C(800, 1)
            integrator.u[1] *= (1 - u2_S) # S
            integrator.u[3] *= (1 - u2_C) # C
        end)
        push!(callbacks, cb)
    end
    
    # Radiotherapy pulse at day 200 + 14*6
    cb_radio = PresetTimeCallback(200.0 + 84.0, (integrator) -> begin
        u1 = get_radiotherapy_effect(60.0)
        integrator.u[3] *= (1 - u1) # C
        integrator.u[4] *= (1 - u1) # Cr
    end)
    push!(callbacks, cb_radio)

    return CallbackSet(callbacks...)
end

function apply_protocol_2(prob)
    # Protocol 2: Protocol 1 + Immunotherapy (10 cycles)
    base_cb = apply_protocol_1(prob)
    
    immuno_cbs = []
    for cycle in 0:9
        start_t = 350.0 + cycle * 20.0
        cb = PresetTimeCallback(start_t, (integrator) -> begin
            # Boost Tc and Th1 populations directly
            dI = 2.0; mtc = 100.0; mth1 = 50.0
            integrator.u[9] += dI * mtc
            integrator.u[7] += dI * mth1
        end)
        push!(immuno_cbs, cb)
    end
    
    return CallbackSet(base_cb, immuno_cbs...)
end

end
