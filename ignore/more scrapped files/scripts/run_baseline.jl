
using DifferentialEquations
using CSV
using DataFrames
include("../src/model.jl")
include("../src/parameters.jl")

using .CancerModel
using .ModelParameters

println("Starting Baseline Simulation...")

p = get_default_parameters()
u0 = get_initial_conditions()
tspan = (0.0, 800.0)

prob = ODEProblem(tumor_ode!, u0, tspan, p)
sol = solve(prob, saveat=1.0)

# Export results for Python plotting
df = DataFrame(sol)
rename!(df, [:time, :S, :Sr, :C, :Cr, :M1, :M2, :Th1, :Th2, :Tc, :Treg, :IL10, :IFNg, :IL2])
CSV.write("baseline_trajectory.csv", df)

println("Simulation complete. Data saved to baseline_trajectory.csv")
