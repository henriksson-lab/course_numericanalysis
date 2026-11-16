# Cross-references from splines

- splines/lec.qmd:24 -> polynomial-interpolation (one polynomial through many points wiggles; About slide)
- splines/lec.qmd:51 -> convergence-cost (logistic OD curve K = 1.2, N0 = 0.05, r = 0.8; rounding noise 1e-16 vs measurement noise; notes)
- splines/lec.qmd:51 -> finite-differences (Numerical differentiation, Part 7, treats the finite-difference side; notes)
- splines/lec.qmd:62 -> polynomial-interpolation (Recap slide: Runge phenomenon on uneven points; notes)
- splines/lec.qmd:134 -> structured-sparse (Thomas algorithm, O(n), no pivoting for diagonally dominant tridiagonal)
- splines/lec.qmd:157 -> structured-sparse (thomas(a, b, c, d) function and its argument convention)
- splines/lec.qmd:239 -> convergence-cost (reading orders off log-log error plots)
- splines/lec.qmd:276 -> convergence-cost (noise amplification of difference quotients, the V; notes)
- splines/lec.qmd:291 -> least-squares-qr (penalised least squares: misfit term; notes)
- splines/lec.qmd:291 -> svd (ridge regularisation, same idea as the roughness penalty; notes)
- splines/lec.qmd:322 -> optimisation-1d (smooth.spline's 1D search over spar; notes)
- splines/lec.qmd:349 -> nonlinear-least-squares (fitting the logistic model directly; notes)
- splines/lec.qmd:415 -> fourier-fft (Next: slide)
- splines/lec.qmd:418 -> finite-differences, fourier-fft (Numerical differentiation, Part 7; Fourier smoothing; notes)
- splines/exercise.qmd:27 -> convergence-cost (OD curve running example)
- splines/exercise.qmd:38 -> structured-sparse (Thomas algorithm, Fact 4)
- splines/exercise.qmd:98 -> structured-sparse (Rust thomas helper comment; same comment in every Rust tab, lines 354, 540, 709, 857, 1042)
- splines/exercise.qmd:187 -> polynomial-interpolation (Newton form, B5(a))
- splines/exercise.qmd:194 -> convergence-cost (V-shaped error plot, B6(c))
- splines/exercise.qmd:206 -> structured-sparse (thomas function, C2)
- splines/exercise.qmd:222 -> polynomial-interpolation (two-point interpolation error formula, Hint 1)
- splines/exercise.qmd:246 -> optimisation-1d (golden-section search, Hint 7)
- splines/exercise.qmd:284 -> convergence-cost (V-shaped error plot, best step and error for central difference, B6 solution)
- splines/exercise.qmd:441 -> structured-sparse (thomas, R code comment in C2 solution; Python at 490)
- splines/exercise.qmd:1392 -> fourier-fft (Where this leads)
- splines/exercise.qmd:1393 -> finite-differences (Where this leads)
- splines/exercise.qmd:1394 -> quadrature (Where this leads: integrating a spline for AUC)
- splines/exercise.qmd:1395 -> nonlinear-least-squares (Where this leads)
- splines/exercise.qmd:1396 -> svd (Where this leads: roughness penalty as regularisation)
- splines/quiz.qmd:18 -> polynomial-interpolation (Q1 hint: Runge)
- splines/quiz.qmd:58 -> structured-sparse (Q4: Thomas algorithm)
- splines/quiz.qmd:106 -> convergence-cost (Q8: V-shaped trade-off)
- splines/quiz.qmd:154 -> nonlinear-least-squares (Q12)

## Running examples defined here

- **Noisy OD reads** (logistic K = 1.2, N0 = 0.05, r = 0.8 from `convergence-cost`, read every 15 min for 10 h, noise N(0, 0.01^2), rounded to 3 decimals; useful for `finite-differences`, `nonlinear-least-squares`, `identifiability`):
  `t <- seq(0, 10, by = 0.25)`;
  `y <- c(0.040, 0.058, 0.076, 0.077, 0.108, 0.127, 0.152, 0.191, 0.200, 0.263, 0.284, 0.327, 0.382, 0.446, 0.502, 0.556, 0.610, 0.672, 0.749, 0.794, 0.839, 0.883, 0.934, 0.958, 1.004, 1.032, 1.076, 1.097, 1.105, 1.111, 1.144, 1.155, 1.163, 1.171, 1.167, 1.182, 1.193, 1.184, 1.177, 1.197, 1.199)`.
  Central differences: rms error of N'(t) 0.0243; smooth.spline (GCV): df 9.416, lambda 1.355e-4 (Python lam 0.1355), rms 0.00666, max rate 0.2388 at 4.01 h (true 0.24 at 3.92 h).
- **Immediate-early gene time course** (TPM; useful for `polynomial-interpolation`, `quadrature`):
  `tg <- c(0, 0.5, 1, 2, 4, 8, 12, 24)`; `eg <- c(1, 40, 32, 6, 1, 0.5, 0.5, 0.5)`. Cubic splines go negative near 3 h (fmm -1.645, natural -1.100); degree-7 polynomial reaches -1.6e5.
- **Gene along pseudotime, 40 cells** (true mean 1 + 2 exp(-((t - 0.6)/0.15)^2), noise sd 0.35):
  `pt <- c(0.013, 0.062, 0.108, 0.126, 0.177, 0.186, 0.202, 0.206, 0.212, 0.266, 0.267, 0.340, 0.372, 0.380, 0.382, 0.384, 0.386, 0.411, 0.482, 0.494, 0.498, 0.573, 0.600, 0.629, 0.652, 0.661, 0.668, 0.687, 0.718, 0.724, 0.770, 0.777, 0.794, 0.827, 0.870, 0.898, 0.908, 0.935, 0.945, 0.992)`;
  `ex <- c(1.32, 1.27, 1.03, 0.30, 1.22, 0.98, 0.95, 0.49, 0.84, 1.16, 1.49, 1.06, 1.33, 1.21, 0.76, 1.11, 1.12, 1.39, 2.46, 2.48, 2.20, 2.85, 3.24, 3.12, 2.53, 2.45, 2.76, 2.70, 2.04, 2.32, 1.69, 1.28, 1.49, 0.81, 1.58, 1.73, 0.90, 0.65, 1.21, 0.95)`.
  GCV: df 9.36, peak at 0.598, height 2.97.
