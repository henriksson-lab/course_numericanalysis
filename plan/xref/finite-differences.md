# Cross-references from finite-differences

- finite-differences/lec.qmd:35 (notes, About this lecture) -> convergence-cost (forward difference of the OD curve, truncation/rounding V, central difference preview)
- finite-differences/lec.qmd:38 (slide "Recap: the V from Convergence and cost") -> convergence-cost (one-slide recap: forward difference, eps/h rounding, V, sqrt(eps) step, difference-ratio check)
- finite-differences/lec.qmd:51 (notes, Recap) -> convergence-cost (students who skipped it should go back)
- finite-differences/lec.qmd:141 (notes, "Any stencil you need") -> lu / Part 3 (the weights solve a small Vandermonde-type linear system)
- finite-differences/lec.qmd:141 (notes, "Any stencil you need") -> diffusion-explicit (same construction builds every stencil in Part 10)
- finite-differences/lec.qmd:153 (slide "The second derivative") -> diffusion-explicit (the (1,-2,1)/h^2 stencil becomes the diffusion term D d2c/dx2; "workhorse of Part 10")
- finite-differences/lec.qmd:158 (notes, The second derivative) -> diffusion-explicit, diffusion-implicit, steady-state-2d, advection-diffusion, reaction-diffusion (Part 10: h^2/12 N'''' is the spatial accuracy of every scheme)
- finite-differences/lec.qmd:183 (notes, All four in R) -> convergence-cost (vectorised over h; difference-ratio check)
- finite-differences/lec.qmd:234 (notes, Failure: second derivative) -> convergence-cost (difference ratios without the exact answer)
- finite-differences/lec.qmd:254 (notes, Richardson extrapolation) -> quadrature (Richardson on the trapezoid rule = Romberg)
- finite-differences/lec.qmd:295 (notes, Measured data) -> error-conditioning (differentiation is ill-conditioned)
- finite-differences/lec.qmd:312 (slide "Smooth first, then differentiate") -> splines (smoothing spline, smooth.spline / make_smoothing_spline, GCV)
- finite-differences/lec.qmd:325 (slide "Smooth first", bullet "fit a growth model ... (Part 9)") -> nonlinear-least-squares (fit a logistic/Gompertz model and differentiate it)
- finite-differences/lec.qmd:328 (notes, Smooth first) -> splines; nonlinear-least-squares
- finite-differences/lec.qmd:364 (notes, Gradients for fitting) -> newton-systems (finite-difference Jacobian, same construction)
- finite-differences/lec.qmd:388 (notes, Failure: optim's default gradient) -> identifiability (Hessian by finite differences for standard errors, Part 9); gradient-methods (BFGS)
- finite-differences/lec.qmd:446 (Summary) -> splines (smooth first)
- finite-differences/lec.qmd:449 (Summary, "Next") -> quadrature (integration averages noise; Romberg)
- finite-differences/lec.qmd:452 (notes, Summary) -> diffusion-explicit (second-derivative stencil returns in Part 10)
- finite-differences/exercise.qmd:26 (Background) -> convergence-cost (running example, logistic OD curve)
- finite-differences/exercise.qmd:271 (Hint 6) -> splines (smoothing spline as library route)
- finite-differences/exercise.qmd:1169 (Where this leads) -> quadrature (Romberg, noise averaging)
- finite-differences/exercise.qmd:1170 (Where this leads) -> splines (smoothing spline for noisy derivatives)
- finite-differences/exercise.qmd:1171 (Where this leads) -> gradient-methods, nonlinear-least-squares, identifiability (gradients/Hessians by finite differences in fitting)
- finite-differences/exercise.qmd:1172 (Where this leads) -> diffusion-explicit (the (1,-2,1)/h^2 stencil of B2 on a spatial grid)

## Note for Part 10 (reuse of the second-derivative stencil)

The second-difference stencil (N(t+h) - 2N(t) + N(t-h))/h^2 with error (h^2/12) N''''
is derived here (lecture slide "The second derivative", exercise B2(a) and C2 via
fd_weights(c(-1,0,1), 2)). Part 10 (*The diffusion (heat) equation by finite
differences* onwards) should refer back to it rather than re-derive it, and can
reuse the rounding warning: its rounding error grows like 4 eps |f| / h^2, best
step ~ eps^(1/4) (lecture slides "Each stencil has its own V" and "Failure: the
second derivative with a tiny step"). The fd_weights() function of exercise C2
(weights from Taylor matching) can build any Part 10 stencil (one-sided boundary
stencils, five-point).

## Running examples defined here

### Simulated plate-reader OD data (noisy logistic growth curve)

The logistic OD curve of convergence-cost (K = 1.2, N0 = 0.05, r = 0.8 /h),
sampled every 10 min from 0 to 10 h (61 points), plus Gaussian noise sd 0.005,
rounded to 3 decimals:

    set.seed(1)
    t <- seq(0, 10, by = 1/6)
    y <- round(od(t) + rnorm(length(t), sd = 0.005), 3)

Values (also hard-coded in exercise C4 for Python/Rust):
0.047 0.058 0.060 0.081 0.084 0.090 0.108 0.123 0.137 0.150 0.177 0.192 0.210
0.226 0.269 0.292 0.322 0.359 0.393 0.428 0.466 0.504 0.540 0.569 0.622 0.659
0.698 0.730 0.772 0.812 0.851 0.876 0.909 0.935 0.955 0.985 1.007 1.029 1.054
1.069 1.079 1.092 1.109 1.120 1.123 1.132 1.145 1.154 1.155 1.166 1.168 1.167
1.175 1.171 1.187 1.192 1.183 1.181 1.191 1.189 1.203

Facts established here (keep consistent):

- True growth rate N'(t) = r N (1 - N/K); max 0.240 OD/h at t = ln(23)/0.8 = 3.92 h.
  Derivatives at t = 3: N = 0.388788, N' = 0.210260, N'' = 0.059212,
  N''' = -0.042270, N'''' = -0.061704.
- Central differences of neighbours: rms error 0.0203 OD/h, range -0.033..0.270.
  Best half-width 40 min (k = 4): rms 0.0053. smooth.spline (GCV): rms 0.0047,
  max 0.2366 at 4.0 h. Local straight-line slope m = 4: rms 0.0041.
- Gaussian log-likelihood with sigma = 0.005, theta = (r, K, N0):
  MLE (BFGS with numDeriv gradient) = (0.7954268, 1.2024571, 0.05076077),
  standard errors (from numDeriv::hessian) = (0.0035, 0.0015, 0.00062).
  optim BFGS with default ndeps = 1e-3 stops at (0.7948541, 1.202606, 0.0508594).
  Candidate reuse: nonlinear-least-squares / identifiability (Part 9), splines.
