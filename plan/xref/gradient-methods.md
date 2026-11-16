# Cross-references from gradient-methods (Gradient descent, Newton and quasi-Newton)

- gradient-methods/lec.qmd:26 (About this lecture) -> newton-systems (Jacobian)
- gradient-methods/lec.qmd:31 (About this lecture) -> intro (growth-curve failure of optim, minutes vs hours)
- gradient-methods/lec.qmd:34 (notes, About) -> optimisation-1d (f' = 0, f'' > 0, golden section, Brent, line search)
- gradient-methods/lec.qmd:57 (notes, biology slide) -> newton (two-site dose-response; solve for dose vs solve for parameters)
- gradient-methods/lec.qmd:83 (slide "What we already have") -> optimisation-1d (minimum conditions, line search)
- gradient-methods/lec.qmd:87 (slide "What we already have") -> newton-systems (Jacobian, Newton step solve(J, -F), FD Jacobian, numDeriv::jacobian)
- gradient-methods/lec.qmd:121 (notes, gradient slide) -> newton-systems (partial derivatives)
- gradient-methods/lec.qmd:146 (notes, gradient of SSE) -> finite-differences (numDeriv Richardson extrapolation)
- gradient-methods/lec.qmd:164 (notes, Hessian) -> intro (Taylor table, quadratic model)
- gradient-methods/lec.qmd:184 (slide "Hessian at the 4PL minimum") -> pivoting-conditioning (condition number)
- gradient-methods/lec.qmd:187 (notes) -> newton-systems (FD Jacobian); identifiability (Hessian -> standard errors)
- gradient-methods/lec.qmd:209 (notes, fixed step) -> euler (stability bound h < 2/lambda; gradient flow)
- gradient-methods/lec.qmd:235 (notes, Armijo) -> optimisation-1d (backtracking line search)
- gradient-methods/lec.qmd:246 (notes, zig-zag) -> pivoting-conditioning (condition number bounds error amplification)
- gradient-methods/lec.qmd:265 (notes, Newton for minimisation) -> optimisation-1d (x - f'/f'')
- gradient-methods/lec.qmd:292 (notes, Newton in R) -> newton (quadratic convergence)
- gradient-methods/lec.qmd:302 (notes, three methods) -> convergence-cost (linear vs quadratic convergence on log axis)
- gradient-methods/lec.qmd:319 (notes, Newton fails) -> newton-systems (Newton finds the unstable middle steady state); nonlinear-least-squares (Levenberg-Marquardt)
- gradient-methods/lec.qmd:324 (slide "Quasi-Newton") -> newton (secant method)
- gradient-methods/lec.qmd:336 (notes) -> newton-systems (Broyden's method)
- gradient-methods/lec.qmd:373 (notes, L-BFGS) -> constrained-lp (constraints)
- gradient-methods/lec.qmd:384 (notes, Nelder-Mead) -> intro (growth-curve plateau)
- gradient-methods/lec.qmd:455 (slide "Back to the growth curve") -> intro (slides "A growth curve, a fit, a green light", "Back to the growth curve: why did optim stop?")
- gradient-methods/lec.qmd:516 (Checklist) -> nonlinear-least-squares (sums of squares)
- gradient-methods/lec.qmd:519 (notes, Checklist) -> identifiability (Hessian, standard errors, near-zero eigenvalues)
- gradient-methods/lec.qmd:530 (Summary, "Next") -> nonlinear-least-squares (Gauss-Newton, Levenberg-Marquardt)
- gradient-methods/exercise.qmd: By hand B6(c) -> newton (secant method)
- gradient-methods/exercise.qmd: In code C3(b) -> newton-systems (FD Jacobian)
- gradient-methods/exercise.qmd: In code C7 -> intro (growth-curve data and failure)
- gradient-methods/exercise.qmd: Hint 4 -> finite-differences (Richardson extrapolation)
- gradient-methods/exercise.qmd: Solution B5(b) -> euler (stability of Euler for y' = -lambda y)
- gradient-methods/exercise.qmd: Where this leads -> nonlinear-least-squares, identifiability, constrained-lp
- gradient-methods/quiz.qmd: Q1 -> newton-systems (Jacobian)

## Running examples defined here

- **4-parameter logistic (4PL) dose-response fit** (agonist reporter assay, simulated
  with set.seed(27), true p = (bottom 5, top 100, EC50 1.5 µM, hill 1.2), noise sd 4,
  rounded to 0.1). Model m(c) = b + (t - b)/(1 + (EC50/c)^h), p = (b, t, EC50, h).
  Doses: rep(10^seq(-2, 2, by = 0.5), each = 2) µM (9 doses x 2 wells).
  y = 12.9, 9.8, 2.9, 0.1, 4.2, 9.7, 17.7, 22.3, 49.7, 42.1, 67.3, 72.6, 97.4, 91.8,
  94.6, 93.3, 92.9, 95.1 (% of max).
  Facts: start p0 = (0, 100, 1, 1), SSE(p0) = 558.9029. Best fit
  p* = (5.393180178, 96.173705843, 1.270980010, 1.239553930), SSE = 314.9239594.
  Hessian eigenvalues at p*: 3708.8, 3618.2, 9.74, 6.74; kappa = 550.
  Steepest descent (Armijo, |g| < 1e-4) 3453 steps; Newton 5; own BFGS 13.
  optim from p0: NM 335 fn; BFGS 40 fn / 11 gr; CG hits maxit (SSE 321.2).
  Dose in nM: kappa 1.9e6; BFGS from EC50 = 1 nM -> flat line (hill -13.8, SSE 25807).
  Newton from (0, 100, 10, 1): Hessian eigenvalues 115.3, 13.1, -42.6, -2422.7.
  Natural reuse: nonlinear-least-squares (Gauss-Newton / LM on the same data),
  identifiability (standard errors from the Hessian; SE and profile of EC50).
