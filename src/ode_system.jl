""" ode_system.jl
Defines the 2D Lotka-Volterra ODE system and related utilities.
The model represents a predator-prey interaction.

System:
    dx/dt = α*x - β*x*y       (prey dynamics)
    dy/dt = -γ*y + δ*x*y      (predator dynamics)

where:  x = prey population
        y = predator population
        α = intrinsic prey growth rate
        β = predation efficiency
        γ = predator mortality rate
        δ = conversion efficiency (predator production per prey consumed)"""

""" lotka_volterra!(du, u, p, t)
In-place ODE function for the 2D Lotka-Volterra system.

Arguments:  du::Vector: Time derivatives (modified in place)
            u::Vector: State variables [x, y] (prey, predator)
            p::NamedTuple: Parameters (α, β, γ, δ)
            t::Float64: Time (not used, included for ODE solver compatibility)"""
function lotka_volterra!(du, u, p, t)
    x, y = u
    α, β, γ, δ = p.α, p.β, p.γ, p.δ
    du[1] = α * x - β * x * y
    du[2] = -γ * y + δ * x * y
end

""" lotka_volterra(u, p, t)
Non-mutating version of the Lotka-Volterra system.
Returns: Vector of time derivatives [dx/dt, dy/dt] """
function lotka_volterra(u, p, t)
    du = similar(u)
    lotka_volterra!(du, u, p, t)
    return du
end

""" has_equilibrium(u, p; tolerance=1e-6)
Check if a state is at an equilibrium point (fixed point).
Arguments:  u::Vector: State variables [x, y]
            p::NamedTuple: Parameters
            tolerance::Float64: Tolerance for checking if derivatives are near zero
Returns:    Bool: true if the state is near an equilibrium """
function has_equilibrium(u, p; tolerance=1e-6)
    du = similar(u)
    lotka_volterra!(du, u, p, 0.0)
    return all(abs.(du) .< tolerance)
end

""" get_equilibrium_points(p)
Calculate the theoretical equilibrium points for the Lotka-Volterra system.
Returns: Tuple: (eq1, eq2) where each is a 2D vector representing an equilibrium point
            eq1 = [0, 0] (extinction)
            eq2 = [γ/δ, α/β] (coexistence) """
function get_equilibrium_points(p)
    α, β, γ, δ = p.α, p.β, p.γ, p.δ
    eq1 = [0.0, 0.0]
    eq2 = [γ / δ, α / β]
    return eq1, eq2
end

""" get_state_summary(u)
Provide a human-readable summary of the current state.
Arguments: u::Vector: State variables [x, y]
Returns: String: Formatted summary """
function get_state_summary(u)
    return "Prey: $(u[1]), Predator: $(u[2])"
end

export lotka_volterra!, lotka_volterra, has_equilibrium, get_equilibrium_points, get_state_summary