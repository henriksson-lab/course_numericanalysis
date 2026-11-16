**Implicit scheme** — a time-stepping method in which the new values appear on both sides of the update, so each step requires solving an equation (here a linear tridiagonal system).
**Backward Euler (for diffusion)** — the second difference evaluated at the new time level: (I - rA) c^{n+1} = c^n + b; first order in time, unconditionally stable, damps fast modes.
**Crank–Nicolson** — the second difference averaged over the old and new time levels (trapezoidal rule in time); second order in time, unconditionally stable, but flips fast modes (growth factor near -1).
**Mesh ratio** — r = D Δt / Δx², the single number through which the time step, grid spacing and diffusion coefficient enter a diffusion scheme.
**Growth (amplification) factor** — the factor by which one time step multiplies a Fourier mode of the grid; |g| ≤ 1 for every mode means the scheme is stable.
**Unconditional stability** — stability for every time step, as for backward Euler and Crank–Nicolson on the diffusion equation.
**Rannacher start** — a few backward Euler steps before switching to Crank–Nicolson, to damp the high-frequency content of sharp initial data.
**Maximum principle** — for pure diffusion, the solution never exceeds the largest (or falls below the smallest) initial or boundary value; a violation signals instability or oscillation.
**Method of lines** — discretising space only, turning a PDE into a large system of ODEs in time that is handed to an ODE solver.
**Banded Jacobian** — a Jacobian whose non-zeros lie near the diagonal (tridiagonal for 1D diffusion); telling a stiff solver about it (`ode.1D`, `jac_sparsity`, `lband`/`uband`) makes each step O(N).
**Exact-solution test** — verifying a solver on a problem with a known solution, such as a sine mode decaying as e^{-D k² t}, and checking the observed order.
**Observed order** — log2 of the ratio of errors (or of successive differences) when a step is halved; should match the design order once in the asymptotic range.
**Asymptotic range** — step sizes small enough that the leading error term dominates and the observed order equals the design order.
**Error floor** — the level at which the error stops decreasing when one step (Δt or Δx) is refined, because the error from the other one dominates.
**Error cancellation** — accidental near-cancellation of time and space errors of opposite sign at one particular step size, making a single run look more accurate than the method is.
**Thermal diffusivity** — D = k/(ρ c_p), the diffusion coefficient of heat; about 1.4×10⁻⁷ m²/s for water, about 1×10⁻⁷ for PDMS and 7×10⁻⁷ for glass.
**Biot number** — Bi = hH/k, the ratio of heat loss at a surface to conduction through a wall of thickness H; sets the steady temperature drop across the wall.

## Install

None beyond intro/install.qmd (uses deSolve in R; numpy, scipy, matplotlib in Python; ode_solvers and std in Rust).
