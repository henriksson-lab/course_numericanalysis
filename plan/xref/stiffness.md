# Cross-references from stiffness

- stiffness/lec.qmd:33 (notes, About this lecture) -> euler (Euler instability on drug clearance); runge-kutta (RK4, adaptive RK45, deSolve::ode / solve_ivp); newton-systems ("Part 2", Newton inside implicit steps); lu ("Part 3", linear solves) [ok]
- stiffness/lec.qmd:57 (notes, A biology problem first) -> runge-kutta (ode45 = adaptive Dormand–Prince RK45) [ok]
- stiffness/lec.qmd:109 (Recap: Euler's stability limit) -> euler (drug clearance C' = -0.8 C, factor 1 - 0.8h, instability slide) [fixed: linked]
- stiffness/lec.qmd:139 (notes, The test equation) -> newton-systems (linearisation with the Jacobian); power-iteration ("Part 5", decoupling into eigenvalue modes) [ok: power-iteration shows the eigen-decomposition argument for L^k x0, close enough for notes]
- stiffness/lec.qmd:233 (Implicit = one nonlinear solve per step) -> newton-systems (Newton for G(y) = 0) [fixed: linked]
- stiffness/lec.qmd:245 (notes, same slide) -> newton-systems ("Part 2"); lu (linear solve per Newton iteration) [ok]
- stiffness/lec.qmd:268 (Backward Euler with Newton in R) -> newton-systems (newton_sys function) [fixed: linked]
- stiffness/lec.qmd:317 (notes, Error against step size) -> convergence-cost (reading slopes / order on log–log plots) [ok]
- stiffness/lec.qmd:334 (notes, Higher order: BDF) -> finite-differences (one-sided second-order difference formula) [ok]
- stiffness/lec.qmd:351 (notes, What one implicit step costs) -> lu (LU cost 2/3 n^3); structured-sparse (banded/sparse Jacobians); newton-systems (finite-difference Jacobian) [ok]
- stiffness/lec.qmd:456 (When implicit methods mislead) -> euler (24-hour oscillator x'' = -omega^2 x, Euler spiral, amplitude factor) [fixed: linked]
- stiffness/lec.qmd:478 (notes, Checks for stiff problems) -> newton-systems (check an analytic Jacobian against a finite-difference one) [ok]
- stiffness/lec.qmd:490 (Summary, "Next") -> phase-plane (Jacobian eigenvalues decide stability of steady states) [fixed: linked]
- stiffness/exercise.qmd:37 (Background, fact 2) -> newton-systems (Jacobian definition) [fixed: linked]
- stiffness/exercise.qmd:141 (B1) -> euler (drug clearance example) [fixed: linked]
- stiffness/exercise.qmd:179 (B7(b)) -> diffusion-implicit (diffusion on a fine grid as a stiff ODE system, method of lines) [fixed: linked]
- stiffness/exercise.qmd:208 (Hint 2) -> newton (one-variable Newton update) [fixed: linked (newton)]
- stiffness/exercise.qmd:212 (Hint 3) -> newton-systems (Jacobian, row = equation) [ok: plain, second mention on page]
- stiffness/exercise.qmd:1300 (Where this leads) -> phase-plane (eigenvalues of the Jacobian and stability of steady states) [fixed: linked]
- stiffness/exercise.qmd:1301 (Where this leads) -> nonlinear-least-squares (fitting stiff ODE models, choice of solver) [fixed: linked]
- stiffness/exercise.qmd:1302 (Where this leads) -> diffusion-implicit (method of lines, one tridiagonal solve per backward Euler step) [fixed: linked]

## Running examples defined here

- **Mass-action enzyme kinetics** E + S <-> C -> E + P with kf = 100 per uM per s, kr = 100 per s, kcat = 1 per s, S0 = 10 uM, E0 = 1 uM, C0 = P0 = 0, simulated over 0–30 s. Fast eigenvalue about -1201 per s (binding, ~1 ms), slow about -0.083 per s (~12 s). Michaelis–Menten: Vmax = 1 uM/s, Km = 1.01 uM. Reference: P(10 s) = 8.002674, P(30 s) = 9.999828. Useful for *Nonlinear least squares and ODE fitting* (fitting kcat, Km) and *Identifiability and uncertainty of fitted parameters* (kf and kr only identifiable through Km when binding is fast).
- **Robertson's problem** (rates 0.04, 3e7, 1e4; y(0) = (1, 0, 0)); reference y2(40) = 9.185535e-6.

## Pass 2

- Overflow of "Try it: the fast eigenvalue sets the step" (741/720 px) fixed: demo embed limited to max-width 1120 px (plots scale with width); re-probed at 720 px.
