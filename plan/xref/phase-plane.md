# Cross-references from phase-plane

- phase-plane/lec.qmd:24 (About this lecture) -> newton-systems (Newton found the three toggle steady states, the middle one most often) [fixed: linked]
- phase-plane/lec.qmd:33 (notes, About this lecture) -> newton-systems (Jacobian, toggle switch); power-iteration, eigen-algorithms (eigenvalues, "Part 5"); euler, runge-kutta (ODE solvers of Part 8) [ok]
- phase-plane/lec.qmd:36-51 (slide "Where we left the toggle switch") -> newton-systems (model, steady states, Newton basin shares 21/57/21 % from fig06-basins.svg) [ok: alpha = 3, n = 2, states and 21/57/21 % match newton-systems]
- phase-plane/lec.qmd:95 (notes, Trajectories with deSolve::ode) -> newton-systems (same F); runge-kutta (adaptive step), stiffness (lsoda switching Adams/BDF) [ok]
- phase-plane/lec.qmd:113 (notes, Nullclines) -> newton-systems (the curves F1 = 0, F2 = 0) [ok]
- phase-plane/lec.qmd:117-128 (slide "Stability in one variable") -> convergence-cost (logistic running example r = 0.8, K = 1.2) [ok]
- phase-plane/lec.qmd:145, 149 (slide "Linearisation in two variables" + notes) -> newton-systems (linear approximation F(x + delta) ≈ F(x) + J delta; entries of J at U on) [fixed: linked]
- phase-plane/lec.qmd:203 (notes, trace and determinant) -> power-iteration (recap: eigenvalues of a 2x2); newton-systems (det J = -0.4186 at the middle state) [ok: det J = -0.4186 in newton-systems exercise B-solution]
- phase-plane/lec.qmd:228 (notes, toggle states classified) -> newton-systems (J); eigen-algorithms (QR algorithm inside eigen()) [ok]
- phase-plane/lec.qmd:272-276 (slide "Eigenvalues of big networks" + notes) -> power-iteration (largest modulus, shifts, Rayleigh quotient); eigen-algorithms (Arnoldi, RSpectra, "largest real part") [fixed: linked]
- phase-plane/lec.qmd:310 (notes, Where Newton goes ≠ where the cell goes) -> newton-systems (basin figure fig06-basins.svg; promise of this contrast) [ok]
- phase-plane/lec.qmd:441 (notes, Failure Re lambda ≈ 0) -> newton-systems (Newton loses quadratic convergence when J is singular at the root, alpha = 2) [ok: newton-systems exercise E4 solution: 39 steps, linear convergence at alpha = 2]
- phase-plane/lec.qmd:485, 488 (Summary, Next + notes) -> gillespie (Stochastic simulation: noise kicks the switch across the separatrix); newton-systems; Part 5 [fixed: linked]
- phase-plane/exercise.qmd:27 (Background) -> newton-systems (toggle switch model and steady states) [fixed: linked]
- phase-plane/exercise.qmd:46 (Fact 8) -> runge-kutta (what the solvers do inside); eigen-algorithms (how eigenvalues are computed) [fixed: linked eigen-algorithms; runge-kutta linked in new Sel'kov sentence]
- phase-plane/exercise.qmd:62, 95, 137 (code comments in Background tabs) -> newton-systems (same Jacobian) [fixed: added "toggle, J from the Background code" comments to C2/C3/C6 tabs that use them]
- phase-plane/exercise.qmd:184 (B4) -> newton-systems (Jacobian formula) [ok]
- phase-plane/exercise.qmd:217 (C3) -> newton-systems (exercise C3: Newton from the same 10 x 10 grid, 23 / 54 / 23 %) [ok: 23/54/23 % matches newton-systems C3 solution]
- phase-plane/exercise.qmd:Where this leads -> gillespie (noise and switching); nonlinear-least-squares (fitting ODE parameters, Part 9); reaction-diffusion (Turing patterns, Part 10) [fixed: linked]
- phase-plane/quiz.qmd:82 (Q6) -> power-iteration, eigen-algorithms (largest modulus vs largest real part; RSpectra / scipy eigs with which = "LR") [ok]

## Running examples defined here

### Biased toggle switch

The toggle switch of newton-systems with a 3 % weaker promoter for V:

    du/dt = 3   / (1 + v^2) - u
    dv/dt = 2.9 / (1 + u^2) - v

Steady states (Newton, 7 digits): U on (2.6555235, 0.3601677), saddle (1.0943960, 1.3195598),
V on (0.4282447, 2.4505795). Saddle eigenvalues 0.176727 and -2.176727. Separatrix
(from backward integration) passes through v = 0.2674, 0.6720, 1.2144, 2.3486, 3.5031 at
u = 0.1, 0.5, 1, 2, 3 (above the diagonal). ODE basins on [0, 4]^2 (300 x 300 grid):
U on 58 %, V on 42 %; Newton basins: 21.8 / 57 / 21.2 %. On the 10 x 10 grid
{0.2, ..., 3.8}: ODE 58 / 0 / 42 (symmetric switch: 45 / 10 / 45, the 10 diagonal starts
converge to the saddle). Bistability is lost for alpha_v between 2.7 and 2.8.

### Facts about the symmetric toggle (alpha = 3, n = 2) established here

- Outer states: eigenvalues -1 ± 2/alpha (pq = 4/alpha^2), i.e. -1/3 and -5/3 at alpha = 3;
  eigenvector for -1/3 is (0.9342, -0.3568), for -5/3 (0.9342, 0.3568).
- Middle state: eigenvalues -1 ± g, g = 2 s^2/(1 + s^2) = 1.1911; unstable direction (1, -1),
  stable direction (1, 1). Separatrix = diagonal u = v.
- Starting exactly on the diagonal, simulations converge to the saddle; from 1e-4 off it takes
  about 50 time units to leave, from 1e-8 about 100.
- At alpha = 2: pitchfork bifurcation at (1, 1), eigenvalues 0 and -2. multiroot from (3, 0.5)
  returns (1.0017321, 0.9982688).

### Sel'kov glycolysis oscillator (limit cycle)

    du/dt = -u + a v + u^2 v      (u = ADP)
    dv/dt =  b - a v - u^2 v      (v = fructose-6-phosphate)

a = 0.08, b = 0.6. Fixed point (b, b/(a + b^2)) = (0.6, 1.363636); J there
[[0.6364, 0.44], [-1.6364, -0.44]], trace 0.1964, det a + b^2 = 0.44, eigenvalues
0.09818 ± 0.65602i (unstable spiral). Limit cycle: u in [0.2453, 1.5176], period 9.6226.
Hopf points (trace = 0) at b = 0.3464 and b = 0.8485 for a = 0.08 (b^4 + (2a - 1) b^2 + a + a^2 = 0).

## Pass 2

- Repressilator: now named with a link to runge-kutta on slide "Oscillations: glycolysis in yeast" and in the exercise Background (six variables, so Sel'kov is used for the phase plane); notes give alpha = 216, period about 8.7.
- Exercise E3 wrongly attributed the biased-switch Newton basins (21.8/57/21.2 %) to newton-systems; now points to this lecture's slide, with 21/57/21 % for the symmetric switch in newton-systems.
