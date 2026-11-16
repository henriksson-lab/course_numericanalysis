# Cross-references from diffusion-explicit

- diffusion-explicit/lec.qmd:33 (notes, About this lecture) -> finite-differences (second-difference stencil), euler (time step), power-iteration / eigen-algorithms (Part 5, eigenvalues), newton-systems (partial derivatives); units 32-35 reuse the chamber
- diffusion-explicit/lec.qmd:44 (notes, Biology slide) -> advection-diffusion (flow-based gradient generators, "Flow plus diffusion")
- diffusion-explicit/lec.qmd:51 (slide "Recap: two tools you already have") -> finite-differences (the (1,-2,1)/h^2 stencil and its h^2/12 error)
- diffusion-explicit/lec.qmd:58 (same slide) -> euler (amplification factor 1 - kh, |1 - kh| <= 1)
- diffusion-explicit/lec.qmd:69 (notes, Recap) -> stiffness (eigenvalue-wise stability for y' = Jy)
- diffusion-explicit/lec.qmd:77 (slide "A function of space and time") -> newton-systems (partial derivatives as slopes of slices)
- diffusion-explicit/lec.qmd:80 (notes) -> newton-systems
- diffusion-explicit/lec.qmd:113 (notes, Fick to diffusion equation) -> euler (rate in - rate out); diffusion-implicit (heat: Drosophila embryo between two temperatures)
- diffusion-explicit/lec.qmd:169 (notes, Step 1 space) -> diffusion-implicit (method of lines); structured-sparse (same tridiagonal matrix as the morphogen in a row of cells)
- diffusion-explicit/lec.qmd:280 (notes, exact solution) -> euler (order check on log-log axes, Part 8)
- diffusion-explicit/lec.qmd:309 (notes, zigzag) -> fourier-fft (Nyquist mode, shortest wave on the grid); euler (drug clearance, 1/k < h < 2/k sign flips)
- diffusion-explicit/lec.qmd:326 (slide "The eigenvalue view") -> power-iteration / eigen-algorithms (Part 5), stiffness (stability for every eigenvalue)
- diffusion-explicit/lec.qmd:329 (notes) -> eigen-algorithms (eigen(); Python numpy.linalg.eigvalsh)
- diffusion-explicit/lec.qmd:338 (slide "Slow modes, fast modes") -> diffusion-implicit (fix for the step limit)
- diffusion-explicit/lec.qmd:341 (notes) -> stiffness (definition: step set by stability, not accuracy); diffusion-implicit (backward Euler, Crank-Nicolson, Thomas)
- diffusion-explicit/lec.qmd:387 (notes, checklist) -> euler (step halving, Part 8)
- diffusion-explicit/lec.qmd:398 (Summary, Next) -> diffusion-implicit
- diffusion-explicit/lec.qmd:401 (notes, Summary) -> diffusion-implicit, structured-sparse (Thomas algorithm), stiffness, steady-state-2d, advection-diffusion, reaction-diffusion
- diffusion-explicit/exercise.qmd:31 (Background item 4) -> finite-differences (second difference and its error)
- diffusion-explicit/exercise.qmd:Background item 5 -> euler (FTCS = Euler on the grid ODEs)
- diffusion-explicit/exercise.qmd:Hint 3 -> euler (|1 - kh| <= 1)
- diffusion-explicit/exercise.qmd:264 (Hint 6) -> quadrature (trapezoid rule, "Numerical integration")
- diffusion-explicit/exercise.qmd:325 (Solution B8(c)) -> diffusion-implicit, stiffness
- diffusion-explicit/exercise.qmd:1066 (Where this leads) -> diffusion-implicit, structured-sparse (Thomas), stiffness
- diffusion-explicit/exercise.qmd:1067 (Where this leads) -> steady-state-2d
- diffusion-explicit/exercise.qmd:1068 (Where this leads) -> advection-diffusion
- diffusion-explicit/exercise.qmd:1069 (Where this leads) -> reaction-diffusion
- diffusion-explicit/quiz.qmd:58 (Q4 answer) -> finite-differences, euler
- diffusion-explicit/quiz.qmd:Q8 answer -> diffusion-implicit ("next unit")
- diffusion-explicit/quiz.qmd:142 (Q11 answer) -> stiffness

## Literature cited

- [G1] Kim, Kim & Jeon (2010), [G3] Frontiers review (2015): lec.qmd slide "Biology: a flow-free chemotaxis chamber" (+ notes), notes of "Sink or wall"; exercise.qmd Background. Only the qualitative claims (diffusion-based vs flow-based devices; no shear; reservoirs equilibrate over time) are attributed to them, stated cautiously. All numbers are standard values, not taken from the papers.

## Running examples defined here

### Flow-free chemotaxis chamber (1D diffusion, Part 10 grid setup)

- Straight channel, length L = 500 µm (variant L = 1000 µm), source reservoir at
  x = 0 with c = 1 (concentration scaled by the source), sink reservoir at x = L
  with c = 0 (Dirichlet); no flow. Initial condition: channel empty, c(x, 0) = 0
  for 0 < x <= L (the source node holds 1 from t = 0).
- Units µm and s. Small molecule D = 1000 µm^2/s (1e-9 m^2/s); protein
  chemokine D = 100 µm^2/s (1e-10 m^2/s); O2 2000 µm^2/s; heat in water
  alpha = 1.4e-7 m^2/s = 1.4e5 µm^2/s.
- Grid: N = 50 intervals, dx = 10 µm, 51 points x = seq(0, L, by = dx),
  49 interior unknowns; r = D dt / dx^2 = 0.4, dt = 0.04 s (small molecule).
  R: `ftcs(c, r, nsteps)` (vectorised, c[i] <- c[i] + r*(c[i-1] - 2c[i] + c[i+1]),
  i = 2:(n-1)); `ftcs_closed` adds a no-flux wall at x = L via ghost point
  (c_N += 2r(c_{N-1} - c_N)). Rust versions in rust/examples/diffusion_explicit_*.rs.
- Exact solution: c = 1 - x/L - sum_n 2/(n pi) sin(n pi x/L) exp(-n^2 pi^2 D t / L^2).
- Facts (keep consistent): at t = 60 s, c(x = 100, 200, 300, 400) = 0.765, 0.543,
  0.343, 0.165 (FTCS r = 0.4). Gradient within 5 % of 1 - x/L at t D/L^2 = 0.2578
  (1611 FTCS steps for any D, L with N = 50, r = 0.4): 64.4 s (small, 500 µm),
  257.8 s (small, 1 mm), 644.4 s = 10.7 min (protein, 500 µm), 2577.6 s = 43 min
  (protein, 1 mm). Within 1 %: t D/L^2 = 0.4206 (105.2 s; 70.1 min protein 1 mm).
  Max error vs exact at t = 30 s, r = 0.4: N = 20: 1.114e-3, 40: 2.786e-4,
  80: 6.969e-5, 160: 1.742e-5, 320: 4.356e-6 (ratio 4.00).
  Second-difference matrix (49 x 49) eigenvalues -3.996053 .. -0.003947;
  times D/dx^2: -39.96 .. -0.0395 per s; Euler limit dt < 0.05005 s (r < 0.5005);
  slowest time constant 25.34 s (= L^2/(pi^2 D) = 25.33 s); stiffness ratio 1013.
  r = 0.52 to t = 30 s: values +-4.1e15.
  Dead-end chamber (wall at x = L): at 64.4 s c(L) = 0.327; c(L) reaches 0.95 at
  328.0 s = 1.312 L^2/D. Both ends closed, bolus c = 1 on 230..270 µm: trapezoid
  total 50 conserved, final uniform 0.1.
- Suggested reuse: diffusion-implicit (same grid, compare BE/CN with these
  numbers, e.g. 5 % time 64.4 s and t = 30 s error table), steady-state-2d
  (chamber extended to 2D), advection-diffusion, reaction-diffusion.
