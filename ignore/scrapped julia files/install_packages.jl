"""
install_packages.jl
Install all required Julia packages for the TIMEdynamics project.
Usage: julia install_packages.jl
"""

import Pkg
packages = ["DifferentialEquations", "CSV","DataFrames", "ArgParse", "Statistics", 
    "LinearAlgebra", "Distributions", "StatsBase", "Random", "ProgressMeter", "Dates"]

println("=" ^ 70)
println("Installing required packages for TIMEdynamics")
println("=" ^ 70)
println()

for pkg in packages
    println("Installing $pkg...")
    Pkg.add(pkg)
    println("✓ $pkg installed successfully\n")
end

println("=" ^ 70)
println("All packages installed successfully!")
println("=" ^ 70)
println()
println("You can now run the scripts:")
println("  julia scripts/simulate.jl --help")
println("  julia scripts/bifurcation_analysis.jl --help")
println("  julia scripts/parameter_estimation.jl --help")
