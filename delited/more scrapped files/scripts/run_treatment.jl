
using DifferentialEquations
using CSV
using DataFrames
include("../src/model.jl")
include("../src/parameters.jl")
include("../src/treatment_logic.jl")

using .CancerModel
using .ModelParameters
using .TreatmentProtocols

println("Running Treatment Protocol 2...")

p = get_default_parameters()
u0 = get_initial_conditions()
tspan = (0.0, 800.0)

prob = ODEProblem(tumor_ode!, u0, tspan, p)
cb = apply_protocol_2(prob)

sol = solve(prob, callback=cb, tstops=[200.0, 350.0], saveat=1.0)

df = DataFrame(sol)
rename!(df, [:time, :S, :Sr, :C, :Cr, :M1, :M2, :Th1, :Th2, :Tc, :Treg, :IL10, :IFNg, :IL2])
CSV.write("treatment_protocol2_trajectory.csv", df)

# Calculate Efficacy Metrics
initial_tumor = sol(200.0)[1:4] |> sum
final_tumor = sol(800.0)[1:4] |> sum
fold_change = final_tumor / initial_tumor
th1_th2_ratio = sol.u[end][7] / sol.u[end][8]

println("--- Results ---")
println("Fold Change: ", fold_change)
println("TH1/TH2 Ratio: ", th1_th2_ratio)
println("Data saved to treatment_protocol2_trajectory.csv")
