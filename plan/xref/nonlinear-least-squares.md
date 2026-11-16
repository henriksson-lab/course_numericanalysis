# Cross-references from nonlinear-least-squares

- nonlinear-least-squares/lec.qmd:33 (notes, About this lecture) -> least-squares-qr (Part 4, linear LS by QR), newton-systems (Part 2, Jacobian), runge-kutta / stiffness (Part 8 ODE solvers), gradient-methods (LM's two limits), identifiability (flip-flop leads there)
- nonlinear-least-squares/lec.qmd:102 (Recap slide) -> least-squares-qr (QR, qr.solve, kappa(A^T A) = kappa(A)^2)
- nonlinear-least-squares/lec.qmd:107 (Recap slide) -> newton-systems (Jacobian, linearisation m(p + d) ~ m + J d)
- nonlinear-least-squares/lec.qmd:134 (Linearise the model) -> newton-systems (finite-difference Jacobian, relative step)
- nonlinear-least-squares/lec.qmd:137 (notes) -> convergence-cost (h ~ sqrt(eps) balances truncation and rounding)
- nonlinear-least-squares/lec.qmd:155 (Gauss-Newton step) -> newton-systems (square J, solve(J, -F))
- nonlinear-least-squares/lec.qmd:158 (notes) -> least-squares-qr (Householder QR, Part 4)
- nonlinear-least-squares/lec.qmd:241 (notes, How fast?) -> gradient-methods (Gauss-Newton as quasi-Newton, Hessian J^T J)
- nonlinear-least-squares/lec.qmd:277 (notes, LM slide) -> gradient-methods (gradient descent, trust regions)
- nonlinear-least-squares/lec.qmd:352 (notes, log-parameters) -> constrained-lp (bounds)
- nonlinear-least-squares/lec.qmd:377 (notes, Fitting an ODE model) -> stiffness (lsoda), runge-kutta (deSolve ode / solve_ivp / Dopri5)
- nonlinear-least-squares/lec.qmd:396 (Inner tolerance vs FD step) -> convergence-cost (V-shaped FD error, h ~ sqrt(noise))
- nonlinear-least-squares/lec.qmd:438 (notes, flip-flop) -> identifiability (structural non-identifiability, same model and data)
- nonlinear-least-squares/lec.qmd:463 (Summary, Next) -> identifiability
- nonlinear-least-squares/lec.qmd:466 (notes, Summary) -> identifiability (standard errors from J^T J, profile likelihood), constrained-lp
- nonlinear-least-squares/exercise.qmd:45 (Background, Fact 7) -> convergence-cost (FD error with noise)
- nonlinear-least-squares/exercise.qmd:209 (Hint 3) -> convergence-cost (V-shaped error curve)
- nonlinear-least-squares/exercise.qmd:213 (Hint 4) -> least-squares-qr (QR least squares in Rust with q_tr_mul)
- nonlinear-least-squares/exercise.qmd:221 (Hint 6) -> runge-kutta (Dopri5 continuous output model)
- nonlinear-least-squares/exercise.qmd:242 (Solution B3) -> least-squares-qr (linear LS with design matrix J)
- nonlinear-least-squares/exercise.qmd:365 (Solution C1) -> error-conditioning (catastrophic cancellation near ka = ke)
- nonlinear-least-squares/exercise.qmd:1507 (Where this leads) -> identifiability
- nonlinear-least-squares/exercise.qmd:1508 (Where this leads) -> gradient-methods
- nonlinear-least-squares/exercise.qmd:1509 (Where this leads) -> constrained-lp
- nonlinear-least-squares/exercise.qmd:1510 (Where this leads) -> stiffness

## Running examples defined here

### Oral one-compartment pharmacokinetic model (for `identifiability` and later)

Gut amount A_g, central amount A_c (mg), concentration C = A_c / V (mg/L):

    A_g' = -ka A_g,   A_c' = ka A_g - ke A_c,   A_g(0) = D = 100 mg, A_c(0) = 0
    C(t) = D ka / (V (ka - ke)) (exp(-ke t) - exp(-ka t))      (Bateman function)

Simulation parameters: **ka = 1 /h, ke = 0.2 /h, V = 20 L** (V is strictly V/F),
additive noise N(0, 0.1^2) mg/L, set.seed(1) in R
(`round(bateman(tt, 1, 0.2, 20) + rnorm(12, 0, 0.1), 2)`), then typed in:

    D  <- 100
    tt <- c(0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8, 10, 12)          # h
    y  <- c(1.02, 1.88, 2.34, 2.98, 3.27, 3.26, 3.17, 2.77, 1.92, 1.23, 1.00, 0.61)  # mg/L
    bateman <- function(t, ka, ke, V) D * ka / (V * (ka - ke)) * (exp(-ke * t) - exp(-ka * t))

Facts established here (keep consistent):

- Least-squares fit: ka = 1.0150721, ke = 0.1906157, V = 20.2162104, RSS = 0.06646744
  (nls from (0.5, 0.1, 10) in 6 iterations; residual SE 0.08594 on 9 df;
  nls standard errors ka 0.06243, ke 0.01040, V 0.71647).
- RSS at the true parameters (1, 0.2, 20): 0.08346328; at start (0.5, 0.1, 10): 93.53426.
- True-parameter tmax = ln(5)/0.8 = 2.0118 h, Cmax = 3.3437 mg/L, AUC = D/(V ke) = 25 mg h/L.
- **Flip-flop** (structural non-identifiability): (ka, ke, V) -> (ke, ka, V ke/ka)
  gives the identical curve. Mirror minimum: ka = 0.1906157, ke = 1.0150721,
  V = 3.7963083, same RSS 0.06646744. True-parameter mirror: (0.2, 1, 4).
- Jacobian (forward differences) at the fit: kappa(J) = 169.93, kappa(J^T J) = 28876;
  in log-parameters kappa(J) = 8.734.
- Gauss-Newton from (3, 0.05, 50) jumps to (-2.2247, 0.1787, 0.0443), RSS 6.7e29;
  nls stops with "singular gradient"; nlsLM / LM converge.
- Multistart (nls.lm, log-params, starts ka in {0.03,...,10}, ke = 1.5 x same, V0 = 10):
  14 -> normal solution, 21 -> mirror, 1 runs off (RSS 64.54 = sum(y^2)).
- ODE fit (deSolve lsoda, rtol = atol = 1e-10, nls.lm with epsfcn = 1e-10): V = 20.2162095,
  agrees to 7 digits. At rtol = atol = 1e-3 the reported RSS is 0.06640724 (below the true min).

### Mini example (by hand): IV bolus, one parameter

C(t) = 10 exp(-k t), t = (1, 2, 4) h, y = (6.0, 3.7, 1.3) mg/L. Minimum k* = 0.5044173,
RSS 0.005236049. Gauss-Newton from k0 = 0.3: 0.4581365, 0.502106, 0.5044131.
