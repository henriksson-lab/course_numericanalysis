# Cross-references from least-squares-qr

- least-squares-qr/lec.qmd:24 -> lu, pivoting-conditioning (Part 3 solved square systems; "Part 3" only, no title)
- least-squares-qr/lec.qmd:133 -> structured-sparse (Cholesky for the symmetric positive definite A^T A)
- least-squares-qr/lec.qmd:165 -> pivoting-conditioning (recap slide: matrix condition number, lose log10 kappa digits, kappa(A, exact = TRUE))
- least-squares-qr/lec.qmd:179 -> svd (stretch factors are singular values; notes)
- least-squares-qr/lec.qmd:200 -> error-conditioning (normal equations unstable: ill-conditioned intermediate step; notes)
- least-squares-qr/lec.qmd:213 -> nonlinear-least-squares (4-parameter logistic calibration curve is nonlinear; notes)
- least-squares-qr/lec.qmd:276 -> convergence-cost (reading slopes off a log-log error plot; notes)
- least-squares-qr/lec.qmd:311 -> lu (back substitution)
- least-squares-qr/lec.qmd:412 -> error-conditioning (catastrophic cancellation in the Householder sign choice; notes)
- least-squares-qr/lec.qmd:469 -> polynomial-interpolation (high-degree polynomials wiggle between points)
- least-squares-qr/lec.qmd:472 -> polynomial-interpolation, splines (notes)
- least-squares-qr/lec.qmd:497 -> svd (Next: slide)
- least-squares-qr/exercise.qmd:274 -> error-conditioning (catastrophic cancellation, B4 solution)
- least-squares-qr/exercise.qmd:1069 -> svd (Where this leads: singular values, rank deficiency, lstsq)
- least-squares-qr/exercise.qmd:1070 -> polynomial-interpolation, splines (Where this leads)
- least-squares-qr/exercise.qmd:1071 -> nonlinear-least-squares (Where this leads: 4PL, Gauss-Newton solves linear LS each step)
- least-squares-qr/exercise.qmd:1072 -> identifiability (Where this leads: standard errors from R)

## Running examples defined here

- **qPCR standard curve** (also usable by `svd`, `nonlinear-least-squares`, `identifiability`):
  `logN <- rep(1:7, each = 3)` (10^1..10^7 copies, triplicates);
  `cq <- c(34.52, 34.69, 34.48, 31.62, 31.37, 31.14, 28.05, 28.10, 28.07, 24.54, 24.90, 24.68, 21.13, 20.81, 21.47, 17.89, 17.90, 18.09, 14.71, 14.67, 14.73)`.
  Fit Cq = 37.981429 - 3.333929 log10 N, efficiency 99.5 %, RSS 0.6266, kappa(A) = 10.40.
- **BCA protein standard curve** (useful for `nonlinear-least-squares` 4PL fit, `splines`, `polynomial-interpolation`):
  `conc <- rep(c(0, 25, 125, 250, 500, 750, 1000, 1500, 2000), each = 2)` (BSA, ug/mL);
  `a562 <- c(0.087, 0.103, 0.150, 0.109, 0.226, 0.229, 0.358, 0.344, 0.601, 0.569, 0.780, 0.789, 0.951, 0.942, 1.300, 1.238, 1.545, 1.532)`.
  Simulated from 0.10 + 2.6 (1 - exp(-conc / 2500)) + N(0, 0.015^2). Quadratic in mg/mL: c = (0.1056877, 0.9981095, -0.1428604); unknown 0.70 -> 657 ug/mL. Cubic in ug/mL: kappa(V) = 5.8e9, solve() on normal equations fails.
