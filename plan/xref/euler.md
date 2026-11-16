# Cross-references from `euler` (Initial value problems and Euler's method)

- euler/lec.qmd:34 (notes, About this lecture) -> runge-kutta (library solvers deSolve::ode / solve_ivp come in the next unit)
- euler/lec.qmd:37-54 (slide "Recap: Euler already failed once") -> intro (slide "3. Instability": drug clearance C' = -0.8 C, h = 3 h, factor -1.4; same code and output)
- euler/lec.qmd:72 (notes, slide "ODEs as models") -> convergence-cost (logistic OD growth curve K = 1.2, N0 = 0.05, r = 0.8)
- euler/lec.qmd:106 (notes, slide "Direction field") -> phase-plane (equilibria and their stability)
- euler/lec.qmd:117 (slide "Euler's method: follow the slope") -> finite-differences (forward difference, rearranged)
- euler/lec.qmd:120 (notes) -> intro (Taylor table row "Euler / Runge-Kutta (Part 8)")
- euler/lec.qmd:120 (notes) -> runge-kutta (same template, smarter slopes)
- euler/lec.qmd:159 (notes, slide "Euler in R") -> runge-kutta (library solvers, f(t, y) convention)
- euler/lec.qmd:222 (notes, slide "Measure the order") -> convergence-cost (order from a log-log error-vs-h plot)
- euler/lec.qmd:249 (slide "Without the exact answer: halve the step") -> runge-kutta (higher order)
- euler/lec.qmd:252 (notes) -> runge-kutta (adaptive step-size control from local error estimates)
- euler/lec.qmd:265 (notes, slide "Instability: back to drug clearance") -> intro (h = 3 example), stiffness (fast rate forces tiny steps)
- euler/lec.qmd:288 (notes, slide "Instability in a nonlinear model") -> bisection-fixedpoint (Ricker model, period doubling / chaos of discrete maps)
- euler/lec.qmd:303 (notes, slide "Systems of ODEs: the SIR epidemic") -> runge-kutta (SIR model reused)
- euler/lec.qmd:355 (notes, slide "Higher-order ODEs") -> phase-plane (nonlinear circadian clock models)
- euler/lec.qmd:376 (notes, slide "Euler on an oscillation") -> runge-kutta (RK4's larger stability region)
- euler/lec.qmd:387 (notes, slide "The spiral") -> phase-plane (phase-plane pictures, periodic orbits)
- euler/lec.qmd:399 (notes, slide "Before you trust an ODE solution") -> runge-kutta (adaptive solvers; tolerances)
- euler/lec.qmd:410 (Summary, "Next") -> runge-kutta (higher order, adaptive steps, deSolve::ode, solve_ivp)
- euler/lec.qmd:413 (notes, Summary) -> runge-kutta, stiffness
- euler/exercise.qmd:46 (Background) -> runge-kutta (deSolve::ode, solve_ivp, ode_solvers crate)
- euler/exercise.qmd:180 (B4(d)) -> stiffness (one fast rate limits the step)
- euler/exercise.qmd:283 (Solution B4(d)) -> stiffness
- euler/exercise.qmd:291 (Solution B6(c)) -> runge-kutta (classical RK4 is order 4)
- euler/exercise.qmd:608 (Solution C2) -> runge-kutta (cost of order 1 as motivation)
- euler/exercise.qmd:1032 (Where this leads) -> runge-kutta
- euler/exercise.qmd:1033 (Where this leads) -> stiffness
- euler/exercise.qmd:1034 (Where this leads) -> phase-plane (equilibria, linearisation, oscillations)
- euler/exercise.qmd:1035 (Where this leads) -> nonlinear-least-squares (fitting ODE parameters)

## Running examples defined here

- **SIR epidemic**: S' = -beta S I / N, I' = beta S I / N - gamma I, R' = gamma I,
  with N = 10 000, beta = 0.3 per day, gamma = 0.1 per day (R0 = beta/gamma = 3,
  10 days infectious), start (S, I, R)(0) = (9990, 10, 0), simulated over
  160 days (time in days). Exact peak (from dI/dS):
  I_max = S0 + I0 - (N/R0)(1 + ln(R0 S0 / N)) = 3007.960706, reached when
  S = N/R0 = 3333, at about day 38.36 (Euler h = 0.001 d). Euler peaks:
  h = 4: 3535.245, h = 2: 3249.225, h = 1: 3127.315 (day 41), h = 0.5: 3066.517,
  h = 0.1: 3019.556, h = 0.01: 3009.117 (order 1). Final S at day 160 with
  h = 0.01 d: 593.88 (R = 9406.11). S + I + R is conserved exactly by any
  Runge–Kutta method, so it is no accuracy check. R code in euler/lec.qmd slide
  "SIR in R".
- **24-hour linear oscillator**: x'' = -omega^2 x, omega = 2 pi / 24 per hour,
  x(0) = 1, x'(0) = 0, exact cos(omega t). As a system y = (x, x'). Euler
  amplitude after 10 days: 2849.977 (h = 1 h), 2.275467 (h = 0.1 h),
  1.085723 (h = 0.01 h); per-step growth sqrt(1 + omega^2 h^2). Natural reuse in
  runge-kutta (RK4 on the same problem) and stiffness (stability region on the
  imaginary axis).
- **Drug clearance** C' = -0.8 C, C0 = 10 mg/L (from intro), and the logistic
  OD curve (from convergence-cost) are used as defined there; logistic Euler
  reference values: error at t = 5 h is -0.2004 (h = 1), -0.01703 (h = 0.1),
  -0.001660 (h = 0.01); C ≈ 0.166 in error ≈ C h; Euler with h = 3 h ends in a
  2-cycle 0.7683 / 1.4317, with h = 2.6 h in 1.0121 / 1.3417.
