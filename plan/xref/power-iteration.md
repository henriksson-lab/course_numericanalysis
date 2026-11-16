# Cross-references from power-iteration

- power-iteration/lec.qmd:33 -> svd (notes, About: eigen() used in passing in the SVD unit)
- power-iteration/lec.qmd:33 -> eigen-algorithms (notes, About: how eigen() works, QR algorithm, Krylov methods)
- power-iteration/lec.qmd:91 -> eigen-algorithms (notes, Recap 2x2: characteristic polynomial not used for n >= 4, root finding of it badly conditioned)
- power-iteration/lec.qmd:195 -> convergence-cost (notes, "How fast?": error vs iteration on a log scale, linear convergence)
- power-iteration/lec.qmd:195 -> bisection-fixedpoint (notes, "How fast?": "the population settles" = "the iteration converges", |g'| rate)
- power-iteration/lec.qmd:295 -> bisection-fixedpoint (notes, slow switching: false convergence when the rate is near 1)
- power-iteration/lec.qmd:321 -> least-squares-qr (notes, Rayleigh quotient as a one-column least-squares problem / normal equation)
- power-iteration/lec.qmd:332 -> eigen-algorithms (notes, Rayleigh slide: leading eigenvector of a correlation matrix = PC1, PCA)
- power-iteration/lec.qmd:353 -> lu (notes, Inverse iteration: never form the inverse, never solve(A) %*% b)
- power-iteration/lec.qmd:377 -> lu (notes, Inverse iteration in R: factor once, solve many times; costs 2/3 n^3 and 2 n^2)
- power-iteration/lec.qmd:424 -> pivoting-conditioning (notes, "Shift exactly at the eigenvalue?": ill-conditioned solve, kappa)
- power-iteration/lec.qmd:447 -> eigen-algorithms (notes, "Other eigenvalues": deflation and the QR algorithm find several eigenvalues)
- power-iteration/lec.qmd:466 -> newton (notes, RQI: cubic vs Newton's quadratic convergence)
- power-iteration/lec.qmd:466 -> eigen-algorithms (notes, RQI: shifts in the QR algorithm)
- power-iteration/lec.qmd:483 -> eigen-algorithms (notes, "Power iteration or eigen()?": Arnoldi/Lanczos, RSpectra, irlba, scipy.sparse.linalg.eigs)
- power-iteration/lec.qmd:508 -> eigen-algorithms (Summary: Next slide)
- power-iteration/exercise.qmd:839 -> eigen-algorithms (C5 solution: symmetric eigen routines faster and more accurate)
- power-iteration/exercise.qmd:1085 -> eigen-algorithms (Where this leads: QR algorithm, Krylov methods, PCA of 100 000 cells)
- power-iteration/exercise.qmd:1086 -> phase-plane (Where this leads: Jacobian eigenvalues decide stability)
- power-iteration/exercise.qmd:1087 -> gillespie (Where this leads: continuous-time Markov chains, stationary distribution)
- power-iteration/exercise.qmd:1088 -> iterative-solvers (Where this leads: Jacobi rate set by the largest eigenvalue of the iteration matrix)
- power-iteration/exercise.qmd:B6(c) solution, C7(c) solution -> power-iteration/lec.qmd slide "Shift exactly at the eigenvalue?" (internal)

## Running examples defined here

- **Fish Leslie matrix** (4 age classes 0-3, post-breeding census, annual step; reuse in
  `eigen-algorithms` and any later population example):
  `L <- matrix(0, 4, 4); L[1, ] <- c(0, 1.5, 4, 5); L[cbind(2:4, 1:3)] <- c(0.3, 0.6, 0.5)`.
  Fertilities F = (0, 1.5, 4, 5), survivals s = (0.3, 0.6, 0.5).
  Eigenvalues 1.1773586, -0.2662888 +/- 0.722404i (modulus 0.7699203), -0.6447810;
  |lambda2/lambda1| = 0.654 (damping ratio 1.53); stable age distribution
  (0.6945, 0.1770, 0.0902, 0.0383). Stocked lake: n0 = (0, 1000, 0, 0); after 30 years
  n = (180810, 46072, 23479, 9971). Power iteration (sum normalisation, residual tol 1e-10)
  from (0, 1, 0, 0) converges in 55 steps.
- **Semelparous variant**: same survivals, first row (0, 0, 0, 12); four eigenvalues of
  modulus 1.08^(1/4) = 1.019427 on a circle; power iteration cycles 0.6, 0.5, 12, 0.3.
- **2 x 2 plant (lecture)**: L = (1, 1; 0.75, 0), eigenvalues 1.5 and -0.5, eigenvectors
  (2, 1) and (2, -3). **2 x 2 plant (exercise)**: L = (0.6, 0.9; 0.8, 0), eigenvalues 1.2 and
  -0.6, eigenvectors (3, 2) and (3, -4). Biennial: (0, 4; 0.25, 0), eigenvalues +/-1;
  early-flowering variant (0.4, 4; 0.25, 0), eigenvalues 1.2198, -0.8198.
- **Cell-state Markov chain** (stem-like S, basal B, luminal L; daily; illustrative, after
  Gupta et al. 2011): `P <- matrix(c(0.80, 0.15, 0.05, 0.05, 0.90, 0.05, 0.05, 0.10, 0.85), 3, byrow = TRUE)`
  (row-stochastic, p_{t+1} = p_t P). Eigenvalues 1, 0.80, 0.75; stationary distribution exactly
  (0.20, 0.55, 0.25); eigenvectors of P^T for 0.8 and 0.75: (0, -1, 1)/sqrt2 and (1, -1, 0)/sqrt2.
  Slow-switching cell line: `Pslow <- diag(3) + 0.01 * (P - diag(3))`, same pi, lambda2 = 0.998;
  power iteration with change < 1e-8 stops after 5406 steps with true error 4.98e-6.
  Natural to reuse in `gillespie` (as a continuous-time chain with rates) and `monte-carlo`.
- **Trait correlation matrix** (body length, body depth, fin length; symmetric):
  `C <- matrix(c(1, 0.4, 0.3, 0.4, 1, 0.2, 0.3, 0.2, 1), 3)`; eigenvalues 1.6074674, 0.8113209,
  0.5812117; leading eigenvector (0.6346, 0.5844, 0.5058).
