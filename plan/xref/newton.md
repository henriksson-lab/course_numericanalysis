# Cross-references from `newton` (Newton and secant methods)

- newton/lec.qmd:33 (notes, About this lecture) -> bisection-fixedpoint (previous unit: bisection is safe but slow)
- newton/lec.qmd:36-47 (slide "Recap: bracketing and fixed points") -> bisection-fixedpoint (bisection, 0.3 digits per step, fixed-point condition |g'(x*)| < 1, stopping rules)
- newton/lec.qmd:47 (notes) -> convergence-cost (Newton for sqrt(2) vs halving, slide "Convergence of iterations: digits per step")
- newton/lec.qmd:86 (notes) -> convergence-cost (checking a derivative with a forward difference)
- newton/lec.qmd:135 (notes) -> convergence-cost (quadratic convergence promised there)
- newton/lec.qmd:163 (slide "Why digits double", last bullet) -> bisection-fixedpoint (fixed-point iteration, |g'(x*)| < 1)
- newton/lec.qmd:200 (notes) -> gradient-methods (maximisation view of the MLE, Newton for optimisation)
- newton/lec.qmd:220 (slide "Newton on the score") -> floating-point (expm1 for 1 - exp(-lambda))
- newton/lec.qmd:223 (notes) -> identifiability (standard errors / observed information)
- newton/lec.qmd:312 (notes, Failure 1) -> intro ("confident wrong answer" of the introduction)
- newton/lec.qmd:339 (notes, Failure 2) -> gradient-methods (line search, Part 9)
- newton/lec.qmd:465 (Summary, "Next") -> newton-systems (Jacobian, linear solve per step)
- newton/lec.qmd:468 (notes) -> newton-systems, gradient-methods (Broyden / quasi-Newton)
- newton/exercise.qmd:177 (B1(d)) -> convergence-cost (square-root iteration x = (x + 2/x)/2)
- newton/exercise.qmd:218 (C1(a)) -> convergence-cost (forward-difference error ~ |f''| h / 2)
- newton/exercise.qmd:1222 (Where this leads) -> newton-systems (Jacobian, Broyden)
- newton/exercise.qmd:1223 (Where this leads) -> finite-differences (finite-difference derivative, choice of h)
- newton/exercise.qmd:1224 (Where this leads) -> gradient-methods (Newton for optimisation, line search, trust region)
- newton/exercise.qmd:1224 (Where this leads) -> identifiability (standard error from -S'(lambda-hat))
- newton/exercise.qmd:1225 (Where this leads) -> optimisation-1d (Brent's minimiser: optimize / minimize_scalar)
- newton/quiz.qmd:195 (Q10) -> bisection-fixedpoint (fixed-point convergence rate |g'(x*)|)

## Running examples defined here

- **Two-site dose-response curve** (EC50 without closed form):
  R(c) = 0.3 c/(0.1 + c) + 0.7 c^2/(25 + c^2), dose c in µM (high-affinity site
  K1 = 0.1 µM, n = 1, weight 0.3; low-affinity site K2 = 5 µM, n = 2, weight 0.7).
  c50 (R = 0.5) = 3.261305265055212. Newton from c0 = 10 converges to the
  meaningless root c = -3.04997. f'(c) = 0.03/(0.1 + c)^2 + 35c/(25 + c^2)^2.
- **Zero-truncated Poisson UMI counts** (MLE without closed form): 100 cells in
  which a gene is detected; counts 1, 2, 3, 4, 5 UMIs in 52, 30, 12, 5, 1 cells
  (sum 173, mean 1.73). Score S(l) = 173/l - 100/(1 - exp(-l)); MLE
  lambda-hat = 1.2184695258301268, SE = 1/sqrt(-S'(lambda-hat)) = 0.1325.
  Newton from lambda0 >= 2.8 diverges. Natural reuse in gradient-methods and
  identifiability (likelihood profile, standard error).
- **Michaelis-Menten in log substrate** u = ln(S/Km): Newton on
  1/(1 + e^-u) - 1/2 is u <- u - sinh(u); 2-cycle at u = ±2.17732 (S = 8.82 Km).
