**Test equation** — the linear ODE $y' = \lambda y$ (with $\lambda$ possibly complex), used to study how a method treats one mode of a linearised system.
**Amplification factor** — the factor $R(z)$, $z = h\lambda$, by which one step of a method multiplies the solution of the test equation; Euler $1 + z$, backward Euler $1/(1 - z)$.
**Stability region** — the set of $z = h\lambda$ in the complex plane where $|R(z)| \le 1$, so errors do not grow from step to step.
**A-stable** — a method whose stability region contains the whole left half-plane, so it is stable for every decaying mode at every step size (backward Euler, BDF2, Radau).
**Stiff problem** — an ODE with widely separated time scales, where the step size of an explicit method is set by stability (the fastest mode) rather than by the accuracy needed for the solution of interest.
**Stiffness ratio** — the ratio of the largest to the smallest $|\text{Re}\,\lambda|$ of the Jacobian over the simulated time; ratios of $10^3$ and more signal stiffness.
**Explicit method** — a method that computes $\mathbf y_{n+1}$ from known values only (Euler, RK4, RK45); its stability region is bounded.
**Implicit method** — a method whose step formula contains the unknown $\mathbf y_{n+1}$ on both sides, so each step requires solving an equation (backward Euler, BDF, Radau).
**Backward Euler** — the implicit step $\mathbf y_{n+1} = \mathbf y_n + h\,\mathbf f(t_{n+1}, \mathbf y_{n+1})$; order 1, A-stable.
**BDF (backward differentiation formula)** — implicit multistep methods that fit a polynomial through past values and match its slope at the new point to $\mathbf f$; BDF1 is backward Euler, orders up to 5 are used in practice.
**Radau method** — an implicit Runge–Kutta method (order 5 in common solvers), A-stable, with several coupled stages per step.
**Newton matrix** — the matrix $I - hJ$ (or $I - \gamma hJ$) whose linear system is solved in each Newton iteration of an implicit step; $J$ is the Jacobian of $\mathbf f$.
**Modified Newton** — Newton iterations that reuse an old Jacobian and its LU factorisation over several iterations or steps, to save cost.
**LSODA** — an ODE solver (ODEPACK) that switches automatically between a non-stiff Adams method and a stiff BDF method; the default of `deSolve::ode`.
**Adams method** — an explicit-type linear multistep method for non-stiff problems, used by LSODA in its non-stiff mode.
**Quasi-steady-state approximation (QSSA)** — replacing a fast variable by its momentary equilibrium value ($C' \approx 0$), which removes the fast time scale; gives Michaelis–Menten kinetics for $E + S \rightleftharpoons C \to E + P$.
**Michaelis–Menten kinetics** — the rate law $v = V_{\max}S/(K_m + S)$ with $V_{\max} = k_{cat}E_0$ and $K_m = (k_r + k_{cat})/k_f$, valid when $E_0 \ll S_0 + K_m$.
**Mass-action kinetics** — reaction rates proportional to the product of the reactant concentrations, e.g. $k_f E S$ for binding.
**Robertson's problem** — a three-species chemical kinetics model (Robertson 1966) with rate constants spanning nine orders of magnitude; the standard test problem for stiff solvers.
**Numerical damping** — the artificial decay of an oscillation caused by a method whose amplification factor on the imaginary axis is below 1 (backward Euler).

## Install

No packages beyond those on `intro/install.qmd` (deSolve, scipy, nalgebra and ode_solvers are already listed).
