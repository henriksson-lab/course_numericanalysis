# Cross-references from diffusion-implicit

- diffusion-implicit/lec.qmd:33 -> diffusion-explicit (chamber, grid, FTCS), structured-sparse (Thomas algorithm), stiffness (backward Euler, stiffness, BDF) — speaker notes, plain text
- diffusion-implicit/lec.qmd:65 -> diffusion-explicit (FTCS limit r <= 1/2), steady-state-2d (stricter limit in 2D) — notes
- diffusion-implicit/lec.qmd:70 -> diffusion-explicit ("Recap: FTCS and its limit" slide, the chamber)
- diffusion-implicit/lec.qmd:79 -> diffusion-explicit (chamber numbers), stiffness (stiffness argument) — notes
- diffusion-implicit/lec.qmd:94 -> stiffness (backward Euler with Newton on the enzyme)
- diffusion-implicit/lec.qmd:97 -> stiffness (y_{n+1} = y_n + h f(y_{n+1}), I - hJ), structured-sparse (diagonal dominance, no pivoting) — notes
- diffusion-implicit/lec.qmd:115 -> structured-sparse (Thomas by hand) — notes
- diffusion-implicit/lec.qmd:131 -> structured-sparse (thomas(a, b, c, d) signature and cost)
- diffusion-implicit/lec.qmd:135 -> structured-sparse (vector convention a[1], c[m] unused) — notes
- diffusion-implicit/lec.qmd:156 -> diffusion-explicit (D, L, N, x, exact() Fourier series) — notes
- diffusion-implicit/lec.qmd:167 -> stiffness (A-stability on the negative real axis) — notes
- diffusion-implicit/lec.qmd:189 -> quadrature (trapezoidal vs rectangle rule; "Numerical integration") — notes
- diffusion-implicit/lec.qmd:235 -> diffusion-explicit (refinement with r fixed shows order 2) — notes
- diffusion-implicit/lec.qmd:271 -> convergence-cost (reading the order off a log-log slope) — notes
- diffusion-implicit/lec.qmd:291 -> quadrature (Richardson extrapolation) — notes
- diffusion-implicit/lec.qmd:300 -> stiffness, runge-kutta, euler (Part 8 ODE solvers, plain "Part 8")
- diffusion-implicit/lec.qmd:305 -> diffusion-explicit (compartment model, eigenvalue slide) — notes
- diffusion-implicit/lec.qmd:331 -> stiffness (ode45 riding its stability limit on the enzyme), diffusion-explicit (ghost points) — notes
- diffusion-implicit/lec.qmd:342 -> fem-outlook (finite volumes / finite elements + method of lines) — notes
- diffusion-implicit/lec.qmd:356 -> reaction-diffusion (decay length sqrt(D/k)), structured-sparse (morphogen in a row of cells, sqrt(D/k) cells) — notes
- diffusion-implicit/lec.qmd:391 -> diffusion-explicit (checklist "Before you trust a diffusion simulation") — notes
- diffusion-implicit/lec.qmd:402 -> steady-state-2d (Next)
- diffusion-implicit/lec.qmd:405 -> steady-state-2d, advection-diffusion ("Flow plus diffusion"), reaction-diffusion (IMEX splitting) — notes
- diffusion-implicit/exercise.qmd:28 -> diffusion-explicit (Background fact 1: the chamber and exact())
- diffusion-implicit/exercise.qmd:31 -> structured-sparse (Background fact 4: thomas)
- diffusion-implicit/exercise.qmd:45 (and every Rust thomas doc comment in the tabs) -> structured-sparse (thomas function)
- diffusion-implicit/exercise.qmd:"Where this leads" -> steady-state-2d, advection-diffusion, reaction-diffusion, iterative-solvers
- (lec and exercise also cite literature [M1], [M2], [H1], [H2] from contentplan.md "Literature sources for Part 10")

## Running examples defined here

- **Chamber, implicit steps** (reuses diffusion-explicit's chamber: L = 500 µm,
  D = 1000 µm²/s, N = 50, dx = 10 µm, source c = 1, sink c = 0, empty at t = 0).
  `be_step(c, r)` / `cn_step(c, r)`: c includes the two boundary values; R/Python/Rust
  versions in diffusion-implicit/exercise.qmd Background and rust/examples/diffusion_implicit_background.rs.
  Reference numbers: 15 BE steps of 4 s (r = 40) give c(100 µm, 60 s) = 0.7583 (exact 0.7650);
  15 CN steps of 4 s give c(10 µm, 60 s) = 1.1775 (oscillation; exact 0.9763).
- **Sine-mode verification**: both ends 0, c(x,0) = sin(pi x/L), error at t = 30 s;
  amplitude exp(-D (pi/L)^2 30) = 0.306. N = 400: BE errors 0.0395, 0.0206, 0.0105, 0.00531
  for dt = 6, 3, 1.5, 0.75 s; CN 1.70e-3, 4.22e-4, 1.04e-4, 2.46e-5; CN floor 1.2e-4 at N = 50.
- **Drosophila embryo between two laminar streams** (after [M1, M2]): 1D along the AP axis,
  L = 500 µm, nodes x_i = i L/N (N = 200, dx = 2.5 µm), thermal diffusivity of water
  D = 1.4e5 µm²/s, no flux at the tips (ghost points), exchange term (T_out - T)/tau with
  T_out = 20 °C (x < L/2), 27 °C (x > L/2), 23.5 at the middle node; tau = R^2/(2.405^2 D)
  = 12.3 ms for R = 100 µm; lambda = sqrt(D tau) = R/2.405 = 41.6 µm; 10-90 % width
  2 lambda ln 5 = 134 µm; steady within ~50 ms. T(200 µm) = 21.0518 °C (N = 200).
  This is our own crude model, not [M2]'s 3D finite-element model. reaction-diffusion may
  reuse it for the [M1] temperature-step experiment on Bicoid.
- **PCR chip wall** (after [H1]): 1D wall, heater at z = 0 steps 65 -> 95 °C, channel side at
  z = H; PDMS k 0.15, rho 970, cp 1460 (D = 0.106 mm²/s); glass k 1.1, rho 2200, cp 700
  (D = 0.714 mm²/s); natural convection h = 10 W/(m² K) to air at 25 °C.
  Time to 94.5 °C (insulated top): PDMS 16.6 / 66.4 s, glass 2.5 / 9.8 s for H = 1 / 2 mm;
  steady drop with heat loss: PDMS 4.4 / 8.2 °C, glass 0.63 / 1.25 °C.

## Reconciliation notes

- The chamber parameters were taken from diffusion-explicit/make-diffusion-explicit-figs.R
  and diffusion-explicit/lec.qmd (its plan/xref file did not exist yet when this unit was
  written). Notation follows that unit: A for the (1,-2,1) matrix, r = D dt/dx^2, N intervals.
- The recap slide quotes diffusion-explicit's eigenvalue range -0.039 to -40 per s.
