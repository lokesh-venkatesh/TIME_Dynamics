Based on the paper, here is the exact computational and mathematical workflow the authors used. You can use this step-by-step breakdown as a blueprint to structure the base scripts in your code repository.

### Step 1: Define the Mathematical Model (ODE System)

**Model Formulation:** Implement a 13-dimensional system of Ordinary Differential Equations (ODEs) to represent the various cellular and cytokine populations.
**Parameters:** Define the 71 parameters that dictate the transition rates, interactions, and feedback loops between the variables.
**Initialization:** Set the initial conditions, assuming a single initial cancer stem cell () while initializing the other tumor sub-populations (, , ) to zero.

### Step 2: Parameter Estimation and Calibration

**MCMC-DRAM Execution:** Write a script to estimate unknown parameters using the MCMC-DRAM (Markov chain Monte Carlo - Delayed Rejection Adaptive Metropolis) algorithm.
**Distribution and Iterations:** Assume a normal prior distribution and run the MCMC for 500,000 iterations to ensure chain convergence.
**Dataset Fitting:** Fit the initial growth kinetics of the model to experimental data from Gastric Cancer cell lines.
**Validation:** Validate the resistant cancer parameters against time-course data from Breast cancer (MCF-7/TAX), Hepatocellular Carcinoma (SK-Hep1/CDDP3), and Colon Cancer (SW-620-L-OHP and LoVo-L-OHP) cell lines.
**Manual Calibration:** Manually adjust remaining parameters within biologically feasible ranges to match multi-source cytometric data (from Gastric, Ovarian, and Osteosarcoma studies) and steady states.

### Step 3: Base Simulation without Interventions

**Temporal Evolution:** Solve the ODEs over an extended period (e.g., 800 days) without therapies.
**Steady States:** Record the dynamics and equilibrium points of the variables during the exponential growth phase and subsequent stabilization.

### Step 4: Global Sensitivity Analysis

**eFAST Algorithm:** Use the extended Fourier Amplitude Sensitivity Test (eFAST) method via a MATLAB toolbox to analyze the model.
**Sampling Setup:** Apply the analysis across all 71 parameters. Configure the algorithm to take 100 samples per search curve and resample 5 times (, ).
**Execution:** Run the resulting 36,000 total model simulations.
**Output:** Calculate Sensitivity Indices (with a significance of ) specifically for the variables driving tumor growth (, , , ).

### Step 5: Treatment Protocol Design and Optimization

**Control Variables:** Augment the ODE system with control variables () representing mathematical perturbations for Radiotherapy, Chemotherapy, and Immunotherapy.
**Simulation Grid:** Simulate 1,000 different treatment combinations by varying dosage, duration, and cycle numbers.
**Timing:** Start therapies at a standard tumor detection time of 200 days into the simulation.
**Evaluation Metrics:** Evaluate the efficacy of each treatment run based on the "Fold Change" (reduction in tumor mass) and the  ratio (optimizing for a threshold ).
**Visualization:** Plot the efficacy of these 1,000 combinations on a 4-dimensional scatter plot.
**Protocol Implementation:** Program complex, multi-stage intervention regimens. For example, the paper's "Protocol 2" workflow is written as: 200 Days of Growth  Chemotherapy (6 cycles)  Radiotherapy  15 days of Relaxation  Chemotherapy (6 cycles)  Immunotherapy (10 cycles).

### Step 6: Mathematical Boundedness and Positivity Checks

**Invariant Region Validation:** Compute limits to formally prove that all solutions initiated in the positive 13-dimensional domain () remain positive and bounded over time.
**Equilibrium Existence:** Validate numerically that positive interior equilibria exist for the equations across the 36,000 parameter sets.

---

To build a repository that logically flows from basic model definition to complex treatment optimization without recalculating heavy steps, we need a strict separation of concerns. All mathematical crunching will be done in Julia (leveraging packages like `DifferentialEquations.jl`, `Turing.jl` for MCMC, and `BifurcationKit.jl`), and all outputs will be serialized to `.csv` or `.json` files. Python will strictly read these files and generate `matplotlib` figures.

Here is the repository structure and the exact sequence of steps to reproduce the paper's results.

### Proposed Repository Structure

```text
cancer_remission_project/
│
├── data/                      # Real-world experimental datasets (Gastric, Breast, Colon, etc.)
├── sim_output/                # Julia-generated CSVs and JSONs (intermediate/final results)
├── figures/                   # Final matplotlib PNG/PDF outputs
│
├── src/
│   ├── julia/                 # Simulation & computation logic
│   │   ├── 00_model_def.jl    # 13-D ODE definitions and baseline parameters
│   │   ├── 01_mcmc_dram.jl    # Parameter estimation 
│   │   ├── 02_base_sim.jl     # Unperturbed tumor growth simulation
│   │   ├── 03_efast_sens.jl   # Global sensitivity analysis
│   │   ├── 04_bifurcation.jl  # Steady-state & bifurcation calculations
│   │   └── 05_protocols.jl    # Treatment combination simulations
│   │
│   └── python/                # Plotting logic
│       ├── plot_mcmc.py       # Plots posterior distributions & chain convergence
│       ├── plot_dynamics.py   # Time-series plots of the 13 variables
│       ├── plot_efast.py      # Bar charts/Tornado plots of sensitivity indices
│       ├── plot_bifurc.py     # 1D/2D Bifurcation diagrams
│       └── plot_efficacy.py   # 4D scatter plots & protocol time-series

```

---

### Logical Workflow of Steps

By following this sequence, you pass parameters and states forward, minimizing redundant computations.

#### Step 1: Parameter Calibration (MCMC-DRAM)

**Why first?** You cannot run accurate base simulations, sensitivity analyses, or treatment protocols without the calibrated parameter set.

* **Julia (`01_mcmc_dram.jl`)**: Load the experimental datasets from `data/`. Define prior distributions for the unknown parameters. Run the MCMC-DRAM algorithm for 500,000 iterations to fit the model's initial growth kinetics to the empirical data. Calculate the maximum a posteriori (MAP) estimates.
* **Action**: Save the final calibrated 71 parameters to `sim_output/calibrated_params.json` and the Markov chains to `sim_output/mcmc_chains.csv`.
* **Python (`plot_mcmc.py`)**: Read the chains. Plot trace plots to prove convergence and corner/pair plots to show the posterior distributions of the fitted parameters.

#### Step 2: Base Simulation & Equilibrium Check

**Why next?** Now that parameters are locked, we simulate the natural progression of the tumor (no treatments) to find its natural steady state.

* **Julia (`02_base_sim.jl`)**: Load `calibrated_params.json`. Set initial conditions (, rest 0). Run the ODE solver over  days.
* **Action**: Save the time-series data for all 13 variables to `sim_output/base_dynamics.csv`. Record the final values (steady-state tumor volume) to `sim_output/base_steady_state.json`.
* **Python (`plot_dynamics.py`)**: Plot the exponential growth phase and the plateau of the cancer stem cells (), resistant cells (), and immune populations over time.

#### Step 3: Global Sensitivity Analysis (eFAST)

**Why next?** Before applying treatments, we must identify which biological mechanisms (parameters) actually drive tumor growth or suppression.

* **Julia (`03_efast_sens.jl`)**: Load `calibrated_params.json`. Set up the eFAST sampling grid (, ). Generate the 36,000 parameter combinations. Loop the ODE solver over all 36,000 sets. Calculate the first-order and total-order sensitivity indices for tumor mass variables.
* **Action**: Save the computed indices and p-values to `sim_output/efast_indices.csv`.
* **Python (`plot_efast.py`)**: Generate bar charts showing the most sensitive parameters (e.g., proliferation rates vs. immune exhaustion rates).

#### Step 4: Bifurcation Analysis

**Why next?** Using the highly sensitive parameters identified in Step 3, we look for mathematical "tipping points" where the tumor switches from uncontrollable growth to remission.

* **Julia (`04_bifurcation.jl`)**: Load base parameters. Using a continuation package (like `BifurcationKit.jl`), pick the top 2 highly sensitive parameters (e.g., cancer proliferation rate vs. CD8+ T-cell activation rate). Compute the equilibrium points as these parameters vary, identifying stable/unstable branches and Hopf/Saddle-node bifurcation points.
* **Action**: Save the branch curves (parameter value vs. steady-state tumor size) to `sim_output/bifurcation_branches.csv`.
* **Python (`plot_bifurc.py`)**: Plot 1D or 2D bifurcation diagrams. Use solid lines for stable equilibria and dashed lines for unstable equilibria, highlighting the threshold needed for tumor eradication.

#### Step 5: High-Throughput Treatment Optimization

**Why next?** Now we know the system's weaknesses (from eFAST) and tipping points (from Bifurcation), we can design therapies.

* **Julia (`05_protocols.jl` - Part 1)**: Activate control variables () in the ODEs. Generate a grid of 1,000 treatment combinations varying the doses and timings of Chemotherapy, Radiotherapy, and Immunotherapy. Start all treatments at  days (using the state from day 200 in Step 2 as the initial condition to avoid recalculating the first 200 days).
* **Action**: For each run, calculate the "Fold Change" (tumor reduction) and the  ratio. Save all 1,000 results to `sim_output/therapy_grid_results.csv`.
* **Python (`plot_efficacy.py` - Part 1)**: Create a 4-dimensional scatter plot (X, Y, Z for therapy doses, color-mapped to Fold Change or  ratio).

#### Step 6: Targeted Protocol Timelines

**Why last?** We select the absolute best regimes from Step 5 and simulate their exact chronological administration to visualize the clinical timeline.

* **Julia (`05_protocols.jl` - Part 2)**: Program the specific multi-stage protocols (e.g., Protocol 2: 6 cycles Chemo  Radio  15-day break  6 cycles Chemo  10 cycles Immuno). Run the solver using `DiscreteCallbacks` or `PresetTimeCallbacks` in Julia to apply the periodic drug pulses.
* **Action**: Save the high-resolution time-series of these specific protocols to `sim_output/protocol_2_timeseries.csv`.
* **Python (`plot_efficacy.py` - Part 2)**: Plot the final tumor mass over time under the specific protocols, clearly marking the windows of active treatment with shaded background regions on the matplotlib figure.

---

# declaring state variables
S, S_R, C, C_R, M_1, M_2, T_H1, T_H2, T_C, T_reg, IL10, IFN_gamma, IL2

# declaring parameters
beta_M2, beta_Tc, beta_Th1CK2, beta_Th1CK3, beta_Th2, beta_Treg, gamma_C, gamma_CR, gamma_M1, gamma_M2, gamma_S, gamma_Tc, gamma_Th1, gamma_Th2, gamma_Treg, delta_C, delta_Ck1, delta_Ck2, delta_Ck3, delta_CR, delta_M1, delta_M2, delta_S, delta_Tc, delta_Th1, delta_Th2, delta_Treg, lambda_M1, lambda_M2, lambda_Tc1, lambda_Tc2, lambda_Tc3, lambda_Tc4, lambda_Th1, lambda_Th2, lambda_Treg2, mu_C1, mu_C2, mu_S, mu_SR, mu_TcS, mu_TcTreg, mu_Th1Ck1, mu_Th1Ck3, mu_TregCk1, C_max, CR_max, k1, k11, k2, k3, k4, k5, k6, k8, k9, ktc1, ktc2, ktc3, ktc4, m_C, m_S, p_1, p_2, r_1, r_2, tck, mu_M1Ck2, mu_M2Ck1, k_7, k_10

# governing equations for the dynamical system

# for the four tumour cell populations
d[S]/dt = (gamma_S*(1-m_S)*(1-p_1-p_2))*S - (delta_S+(p_2*gamma_S)+gamma_S*m_S*p_1/2)*S - (mu_S*S*IFN_gamma)/(k1 + IFN_gamma) - (tck*S*T_C)/(ktc1+T_C)
d[S_R]/dt = (gamma_S*(1-p_1-p_2) - (delta_S+(p_2*gamma_S)))*S_R + m_S*gamma_S*(1-p_1/2-p_2)*S - (mu_SR*S_R*IFN_gamma)/(k2 + IFN_gamma) - (tck*S_R*T_C)/(ktc2+T_C)
d[C]/dt = gamma_C*(1-m_C)*log((C_max)/(C+r_1))*C + gamma_S*(p_1+p_2)*S - delta_C*C - m_C*gamma_C*C + (mu_C1*C*IL10)/(IL10+k3) - (mu_C2*C*IFN_gamma)/(IFN_gamma+k4) - (tck*C*T_C)/(ktc3+T_C)
d[C_R]/dt = gamma_C*C_R*log((CR_max)/(C_R+r_2)) + gamma_S*S_R*(p_1+p_2) + m_C*gamma_C*C - delta_CR*C_R + (mu_C1*C_R*IL10)/(IL10+k5) - (mu_C2*C_R*IFN_gamma)/(IFN_gamma+k6) - (tck*C_R*T_C)/(ktc4+T_C)
- NOTE: 0.5*K_tumor = C_max = CR_max

# for the effector populations:
d[M_1]/dt = gamma_M1*M_1*((C+C_R)/(M_1+lambda_M1)) - delta_M1*M_1 + ((mu_M1Ck2*M1*IFN_gamma)/(IFN_gamma+k_7))
d[M_2]/dt = gamma_M2*M_2*((C+C_R)/(M_2+lambda_M2)) - delta_M2*M_2 + ((mu_M2Ck1*M2*IL10)/(IL10+k_10))
d[T_H1]/dt = gamma_Th1*((T_H1*M_1)/(lambda_Th1+T_H1)) - delta_Th1*T_H1 - (mu_Th1Ck1*IL10*T_H1)/(IL10+k8) + (mu_Th1Ck3*IL2*T_H1)/(IL2+k9)
d[T_H2]/dt = gamma_Th2*((T_H2*M_2)/(lambda_Th2+T_H2)) - delta_Th2*T_H2
d[T_C]/dt = gamma_Tc*T_C*((C+C_R)/(T_C+lambda_Tc1)) + gamma_Tc*((T_C*T_H1)/(T_C+lambda_Tc4)) - mu_TcS*T_C*((S+S_R)/(T_C+lambda_Tc2)) - delta_Tc*T_C - mu_TcTreg*T_C*((T_reg)/(lambda_Tc3+T_reg))

# for the inhibitor populations:
d[T_reg]/dt = gamma_Treg*((T_reg*M_2)/(T_reg+lambda_Treg2)) - delta_Treg*T_reg + mu_TregCK1*((IL10*T_reg)/(T_reg+k11))
d[IL10]/dt = beta_M2*M2 - delta_Ck1*IL10 + beta_Treg*T_reg + beta_Th2*T_H2
d[IFN_gamma]/dt = beta_Th1CK2*T_H1 + beta_Tc*T_C - delta_Ck2*IFN_gamma
d[IL2]/dt = beta_Th1CK3*T_H1 - delta_Ck3*IL2

# parameter values:
DEFAULT_PARAMETERS = (
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