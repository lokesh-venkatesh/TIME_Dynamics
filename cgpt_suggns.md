This is a great project — and at 13 coupled ODEs, you’re absolutely right to think carefully about structure before layering in phase planes, eFAST, and MCMC-DRAM.

What you’re building is really a **modeling platform**, not just a simulation script. So the key idea is:

> Separate the *model definition* from the *analysis workflows*.

If you do this right, everything (phase plane, eFAST, MCMC-DRAM) becomes a thin layer built on the same core engine.

---

# 🔷 1. First Principle: Modular Architecture

You want 4 logical layers:

```
1. Core Model Layer
2. Simulation Layer
3. Analysis Layer
4. Experiment / Run Layer
```

Each layer depends only downward — never sideways.

---

# 🔷 2. Recommended Repository Structure

Here’s a clean, scalable layout:

```
tumor_immune_model/
│
├── src/
│   ├── model/
│   │   ├── parameters.py
│   │   ├── initial_conditions.py
│   │   ├── ode_system.py
│   │   └── model_wrapper.py
│   │
│   ├── simulation/
│   │   ├── solver.py
│   │   └── simulation_runner.py
│   │
│   ├── analysis/
│   │   ├── phase_plane.py
│   │   ├── bifurcation.py
│   │   ├── efast.py
│   │   └── mcmc_dram.py
│   │
│   ├── utils/
│   │   ├── plotting.py
│   │   ├── sampling.py
│   │   └── io.py
│
├── data/
│   ├── experimental/
│   ├── processed/
│
├── experiments/
│   ├── run_baseline.py
│   ├── run_efast.py
│   ├── run_mcmc.py
│
├── results/
│
└── README.md
```

Now I’ll explain *why* this structure matters.

---

# 🔷 3. The Core Model Layer (Most Important Part)

This is sacred. Everything else calls this.

## `ode_system.py`

Contains **only the equations**:

```python
def tumor_immune_ode(t, y, params):
    # unpack variables
    T, I, ... = y
    
    # unpack parameters
    r = params["r"]
    k = params["k"]
    ...
    
    dTdt = ...
    dIdt = ...
    
    return [dTdt, dIdt, ...]
```

Nothing else.

No plotting.
No data loading.
No sensitivity code.

Just math.

---

## `parameters.py`

Define:

```python
def get_default_parameters():
    return {
        "r": 0.5,
        "k": 1.2,
        ...
    }
```

Also define:

```python
def get_parameter_bounds():
    return {
        "r": (0.1, 1.0),
        "k": (0.5, 2.0),
    }
```

eFAST and MCMC will use these bounds.

---

## `model_wrapper.py`

This gives a unified interface:

```python
def simulate(params, t_span, y0):
    ...
    return solution
```

This wrapper is what:

* eFAST calls
* MCMC calls
* phase plane analysis calls

This abstraction is critical.

---

# 🔷 4. Simulation Layer

This handles numerical integration only.

## `solver.py`

Wrap SciPy (or your solver of choice):

```python
def solve_ode(ode_func, params, t_span, y0):
    sol = solve_ivp(...)
    return sol
```

This lets you swap solvers later without touching analysis code.

---

# 🔷 5. Phase Plane & Bifurcation

These are deterministic structural analyses.

## Phase Plane

* Choose 2 variables
* Hold others fixed
* Compute vector field
* Plot nullclines
* Plot trajectories

This module should import:

```python
from model.ode_system import tumor_immune_ode
```

NOT the MCMC code. Keep separation clean.

---

## Bifurcation

You’ll:

* Choose a bifurcation parameter
* Loop over values
* Find steady states (root finding)
* Track stability (Jacobian eigenvalues)

This deserves its own file because it’s conceptually different from time simulation.

---

# 🔷 6. eFAST Module Structure

`analysis/efast.py`

This module should do:

### Step 1: Generate sinusoidal samples

### Step 2: Call model_wrapper.simulate()

### Step 3: Extract scalar output metric

(e.g., tumor size at final time)

### Step 4: Perform FFT

### Step 5: Compute sensitivity indices

Very important design choice:

### You must define a single "model output function"

Example:

```python
def model_output(params):
    sol = simulate(params, t_span, y0)
    return sol.y[0, -1]  # tumor at final time
```

eFAST works on scalar outputs.

If you want multiple outputs:

* Run separate eFAST analyses per output

---

# 🔷 7. MCMC-DRAM Module Structure

`analysis/mcmc_dram.py`

This module should contain:

### 1️⃣ Likelihood function

```python
def log_likelihood(params, data):
    sol = simulate(params, t_span, y0)
    model_prediction = ...
    return -0.5 * np.sum((data - model_prediction)**2 / sigma**2)
```

### 2️⃣ Prior function

```python
def log_prior(params):
    ...
```

### 3️⃣ Posterior

```python
def log_posterior(params, data):
    return log_prior(params) + log_likelihood(params, data)
```

### 4️⃣ DRAM sampler implementation

This should NOT contain ODE equations directly — only calls to `simulate()`.

---

# 🔷 8. Experiments Folder (Where You Actually Work)

This is your sandbox.

## `run_baseline.py`

* Load default parameters
* Run simulation
* Plot time courses

## `run_efast.py`

* Load parameter bounds
* Call eFAST module
* Save sensitivity indices

## `run_mcmc.py`

* Load data
* Define unknown parameters
* Run DRAM
* Save posterior samples

You work mostly in the `experiments/` folder.

Everything else stays stable.

---

# 🔷 9. Workflow Strategy (How You Actually Proceed)

Here’s the correct order:

---

## Phase 1: Stabilize Core Model

1. Implement ODEs
2. Verify solver works
3. Check units
4. Plot trajectories
5. Validate qualitative behavior

Only move forward once stable.

---

## Phase 2: Phase Plane & Bifurcation

* Choose reduced subsystems
* Analyze equilibrium structure
* Identify critical parameters

This helps you understand:

* What parameters are structurally important

---

## Phase 3: eFAST

Use wide parameter bounds.

Goal:

* Rank parameters by influence
* Identify insensitive ones

Outcome:

* Fix low-impact parameters
* Reduce dimension before MCMC

This is important — it dramatically stabilizes MCMC.

---

## Phase 4: MCMC-DRAM

Now:

* Only estimate sensitive parameters
* Fix others
* Use posterior diagnostics

  * trace plots
  * autocorrelation
  * effective sample size

---

# 🔷 10. Critical Design Decisions You Must Make Early

### A. Define outputs clearly

For sensitivity & MCMC:

* Tumor size at final time?
* Time to peak?
* Area under curve?

Be explicit.

---

### B. Parameter vectorization

Internally, MCMC works best with arrays:

```python
theta = np.array([...])
```

But your model prefers dicts.

Solution:

Write conversion utilities:

```python
def vector_to_params(theta):
    ...

def params_to_vector(params):
    ...
```

This avoids chaos later.

---

### C. Reproducibility

* Set random seeds
* Save parameter sets
* Save MCMC chains
* Log experiment metadata

---

# 🔷 11. Common Mistakes to Avoid

❌ Mixing plotting inside ODE function
❌ Hard-coding parameters inside equations
❌ Writing eFAST that directly edits ODE file
❌ Letting MCMC modify global variables
❌ Not isolating simulation wrapper

---

# 🔷 12. High-Level Project Philosophy

Think of your project like this:

* The ODE model = physics engine
* eFAST = sensitivity microscope
* MCMC-DRAM = parameter inference engine
* Phase plane = geometry analyzer

All tools operate on the same engine.

---

# 🔷 13. If I Were Structuring This Personally

I would:

1. Build simulation engine first
2. Make it callable with arbitrary parameter dict
3. Write a clean `model_output()` function
4. Implement eFAST next
5. Use eFAST results to prune parameters
6. Then implement DRAM

This order prevents:

* Identifiability disasters
* Wasted MCMC time
* Massive posterior correlations

---

# 🔷 Final Big Picture

Your instinct was correct — but the key upgrade is:

> Don’t write “secondary” and “tertiary” scripts.
> Write independent analysis modules that all call the same simulation engine.

That’s what makes the project clean, scalable, and publishable.

---

If you’d like, next we can:

* Design your parameter vector mapping scheme
* Sketch pseudocode for eFAST tailored to 13 ODEs
* Or sketch DRAM pseudocode customized for stiff biological systems
* Or discuss how bifurcation results should inform MCMC priors

Where do you want to go next?
