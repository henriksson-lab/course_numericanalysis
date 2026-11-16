# Cross-references from runge-kutta

- runge-kutta/lec.qmd:38 (slide "Recap: Euler's method") -> euler (Euler's method, order 1, logistic errors) [fixed: euler measured logistic error at t = 5 h, here t = 10 h; notes now say so]
- runge-kutta/lec.qmd:slide "Better slopes: midpoint and Heun" (notes) -> quadrature (midpoint and trapezoid rules) [ok]
- runge-kutta/lec.qmd:slide "The classical Runge–Kutta method (RK4)" (notes) -> quadrature (Simpson's rule) [ok]
- runge-kutta/lec.qmd:slide "Reading the order off a log–log plot" (notes) -> convergence-cost (asymptotic regime, order from log–log slope) [ok]
- runge-kutta/lec.qmd:slide "Predator–prey: watch the conserved quantity" (notes) -> euler (Euler spirals out on the 24-hour oscillator) [ok]
- runge-kutta/lec.qmd:slide "The repressilator" (notes) -> phase-plane (limit cycle, unstable steady state) [fixed: phase-plane uses toggle switch + Sel'kov, no repressilator; notes reworded]
- runge-kutta/lec.qmd:slide "deSolve: rk4 versus ode45" (notes) -> euler (exact SIR peak formula I_max) [ok]
- runge-kutta/lec.qmd:400 (slide "Events: stop at the epidemic peak") -> stiffness (lsoda explained there) [ok]
- runge-kutta/lec.qmd:410 (slide "How events are found, and missed") -> bisection-fixedpoint, newton (root finding on the dense output) [ok]
- runge-kutta/lec.qmd:425 (slide "Checklist") -> stiffness (huge number of steps = stiffness) [ok]
- runge-kutta/lec.qmd:439 (slide "Summary", Next) -> stiffness [ok]
- runge-kutta/exercise.qmd:42 (Background, running examples) -> euler (SIR epidemic parameters) [ok]
- runge-kutta/exercise.qmd:Background (running examples) -> convergence-cost (logistic OD curve), euler/intro (drug clearance C' = -0.8 C, C0 = 10) [ok]
- runge-kutta/exercise.qmd:Where this leads -> stiffness (stiff problems, LSODA) [ok]
- runge-kutta/exercise.qmd:Where this leads -> phase-plane (predator–prey orbits, repressilator limit cycle) [fixed: same; reworded to say phase-plane studies limit cycles on toggle/Sel'kov]
- runge-kutta/exercise.qmd:Where this leads -> nonlinear-least-squares (ODE fitting and solver tolerances) [ok]

## Running examples used here (defined elsewhere)

- SIR epidemic exactly as in euler (N = 10 000, beta = 0.3/day, gamma = 0.1/day,
  (S, I, R)(0) = (9990, 10, 0), I_max = 3007.961). Values found here: peak
  time 38.35531 days (event, lsoda rtol 1e-10), S = 3333.333, I = 3007.9607,
  R = 3658.706 at the peak; I(60) = 998.2012187.
- Logistic OD curve (K = 1.2, N0 = 0.05, r = 0.8/h): exact N(1) = 0.1058708,
  N(10) = 1.1908121. One step h = 1 from 0: Euler 0.0883333, midpoint 0.1021440,
  Heun 0.1018991, RK4 0.1057646.
- Drug clearance C' = -0.8 C, C0 = 10, h = 0.5: Euler 6, midpoint/Heun 6.8,
  RK4 6.704, exact 6.70320.

## Running examples defined here

- **Lotka–Volterra predator–prey** (hares x, lynx y): x' = a x - b x y,
  y' = c x y - d y with a = 1, b = 0.1, c = 0.075, d = 1.5 per year,
  (x, y)(0) = (10, 5). Conserved V = c x - d ln x + b y - a ln y; V(10, 5) =
  -3.81332, equilibrium (20, 10) with V = -4.29618. Period about 5.5 years;
  lynx range 3.09–23.28, prey range about 8–40.6. Drift of V over 50 years:
  Euler h = 0.01: 0.5107; Heun h = 0.01: 6.27e-5; RK4 h = 0.1: 5.69e-6.
  Lynx first reach 20 at t = 3.01715 years.
- **Repressilator** (Elowitz & Leibler 2000, dimensionless):
  m_i' = -m_i + alpha / (1 + p_j^n) + alpha0, p_i' = -beta (p_i - m_i),
  (i, j) = (1,3), (2,1), (3,2); alpha = 216, alpha0 = 0.216, beta = 5, n = 2,
  start m = (1, 0, 0), p = (0, 0, 0). Limit cycle, proteins 1–60, period about
  8.7 time units. Reference p1(50) = 47.335989 (tight tolerance). Fixed-step
  RK4: h = 1 blows up (3e58), h = 0.2 gives 58.785, h = 0.01 gives 47.338477.
  Natural reuse in phase-plane (limit cycle) and gillespie (stochastic version).

## Pass 2

- All italic unit titles turned into links.
- [fixed] Python C2, C3, C5, C6 solution tabs now self-contained (steppers / rk_solve / adapt_heun repeated); rerun with python3 -I, outputs unchanged. R tabs keep their one-line "uses ... from C1/C4" pointers (allowed by the brief).
- SIR parameters agree with euler (N 10 000, beta 0.3, gamma 0.1, (9990, 10, 0)).
