
# 13-D Cancer Dynamics Simulation (Julia Implementation)

This repository provides a modular implementation of the 13-dimensional ODE model for tumor-immune interactions described in Ganguli & Sarkar (2018).

## Structure
- `src/`: Core logic (ODE definition, Parameters, Sensitivity, Inference).
- `scripts/`: Entry points for simulation and analysis.
- `plotting/`: Python scripts for visualization.

## Instructions
1. **Julia Dependencies**: Ensure you have installed `DifferentialEquations`, `GlobalSensitivity`, `AdaptiveMCMC`, `DataFrames`, and `CSV`.
2. **Run Baseline**: `julia scripts/run_baseline.jl` to generate `baseline_trajectory.csv`.
3. **Run Treatment**: `julia scripts/run_treatment.jl` to simulate complex therapeutic protocols.
4. **Sensitivity**: Use `src/sensitivity_analysis.jl` functions to identify key tumor drivers.
5. **Inference**: Use `src/mcmc_inference.jl` with your experimental CSV data for parameter fitting.
6. **Visualize**: Run `python plotting/visualize.py` to generate plots from the generated CSV files.

## Model Summary
- **S**: Cancer Stem Cells
- **Sr**: Resistant Cancer Stem Cells
- **C/Cr**: Differentiated Cancer Cells (Sensitive/Resistant)
- **Immune**: M1, M2, TH1, TH2, Tc, Treg
- **Cytokines**: IL10, IFNg, IL2
