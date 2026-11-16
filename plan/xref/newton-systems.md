# Cross-references from newton-systems

- newton-systems/lec.qmd:33 (notes, About this lecture) -> newton (one-variable Newton method, the "previous unit") [ok]
- newton-systems/lec.qmd:33 (notes, About this lecture) -> convergence-cost (quadratic convergence, the sqrt(eps) step) [ok]
- newton-systems/lec.qmd:33 (notes, About this lecture) -> phase-plane (Jacobian for stability) [ok]
- newton-systems/lec.qmd:33 (notes, About this lecture) -> gradient-methods (Jacobian/Hessian in optimisation, "Part 9") [ok]
- newton-systems/lec.qmd:48 (notes, toggle switch slide) -> phase-plane (reuses the toggle-switch model with alpha = 3, n = 2) [ok]
- newton-systems/lec.qmd:72 (notes, steady states slide) -> phase-plane (nullclines; instability of the middle state) [ok]
- newton-systems/lec.qmd:104 (Recap: Newton in one variable) -> newton (sqrt(2) example, tangent-line derivation) [fixed: recap table now uses newton's two-site dose–response EC50 run from c0 = 1 (errors 2.26, 0.36, 6.5e-3, 1.6e-6, 9.1e-14), linked]
- newton-systems/lec.qmd:116 (notes, Recap) -> newton (derivation, failure modes, secant method) [ok]
- newton-systems/lec.qmd:149 (notes, partial derivatives as slopes) -> convergence-cost (forward difference) [ok]
- newton-systems/lec.qmd:247 (notes, one step by hand) -> lu (why never compute the inverse) [ok]
- newton-systems/lec.qmd:270 (notes, Newton for systems in R) -> pivoting-conditioning (norms in general) [ok]
- newton-systems/lec.qmd:297 (notes, Still quadratic) -> convergence-cost (reading linear vs quadratic convergence on a log axis) [ok]
- newton-systems/lec.qmd:329 (notes, basins slide) -> phase-plane (stability from eigenvalues of the Jacobian) [ok]
- newton-systems/lec.qmd:361 (No formula for J? Finite differences) -> convergence-cost (h = sqrt(eps) balances truncation and rounding) [ok, linked]
- newton-systems/lec.qmd:364 (notes, finite differences) -> convergence-cost (forward difference); gradient-methods (quasi-Newton / Broyden, "Part 9"); finite-differences (central differences) [ok]
- newton-systems/lec.qmd:395 (notes, FD Jacobian in R) -> convergence-cost (the V-shaped error vs h plot) [ok]
- newton-systems/lec.qmd:412 (notes, cost of a step) -> lu (LU factorisation, cost, not the inverse); structured-sparse (sparse Jacobians of gene networks) [ok]
- newton-systems/lec.qmd:464 (Failure 2: checks) -> phase-plane (eigenvalues of J at a root decide stability) [ok, linked]
- newton-systems/lec.qmd:467 (notes, Failure 2) -> phase-plane (the middle state is a saddle) [ok]
- newton-systems/lec.qmd:493 (Summary, "Next") -> lu (what solve does inside) [ok, linked]
- newton-systems/lec.qmd:496 (notes, Summary) -> phase-plane (stability of the steady states); gradient-methods (Hessian, Newton for optimisation); lu (Part 3 opens solve) [ok]
- newton-systems/exercise.qmd:279 (Solution B8(b), "Part 3") -> lu / structured-sparse (cost of dense vs sparse solves) [ok, linked]
- newton-systems/exercise.qmd:794 (Solution C4, "Part 9") -> gradient-methods (quasi-Newton reuses an approximate Jacobian) [ok, linked]
- newton-systems/exercise.qmd:1281 (Where this leads) -> lu (what solve(J, -F) does, 2/3 n^3 cost); pivoting-conditioning (norms, nearly singular J amplifies errors) [ok, linked]
- newton-systems/exercise.qmd:1282 (Where this leads) -> structured-sparse (sparse Jacobians of gene networks) [fixed: structured-sparse has no Jacobians; reworded to "sparse storage and solvers; a gene network's Jacobian is sparse"]
- newton-systems/exercise.qmd:1283 (Where this leads) -> finite-differences (central differences, Richardson extrapolation as in numDeriv) [ok, linked]
- newton-systems/exercise.qmd:1284 (Where this leads) -> phase-plane (same toggle switch; eigenvalues of J decide stability) [ok, linked]
- newton-systems/exercise.qmd:1285 (Where this leads) -> gradient-methods (Newton on the gradient, Hessian; quasi-Newton) [ok, linked]
- newton-systems/quiz.qmd:106 (Q8 answer) -> phase-plane (stability via eigenvalues of J) [ok, linked]

## Running examples defined here

### Genetic toggle switch (symmetric, dimensionless)

Two genes U and V whose proteins repress each other (Gardner, Cantor & Collins,
Nature 2000; LacI/TetR in E. coli). Concentrations u, v in units of the
repression threshold; time in units of the protein lifetime:

    du/dt = alpha / (1 + v^n) - u
    dv/dt = alpha / (1 + u^n) - v

Default parameters: **alpha = 3** (promoter strength), **n = 2** (Hill coefficient).

R code used in this unit (Python/Rust equivalents in the exercise):

    alpha <- 3; n <- 2
    F <- function(x) { u <- x[1]; v <- x[2]
      c(alpha / (1 + v^n) - u, alpha / (1 + u^n) - v) }
    J <- function(x) { u <- x[1]; v <- x[2]
      matrix(c(-1, -alpha * n * v^(n - 1) / (1 + v^n)^2,
               -alpha * n * u^(n - 1) / (1 + u^n)^2, -1), nrow = 2, byrow = TRUE) }

Steady states for alpha = 3, n = 2 (exact):

- "U on":  (u, v) = ((3 + sqrt 5)/2, (3 - sqrt 5)/2) = (2.618033988749895, 0.381966011250105)
- "V on":  (0.381966011250105, 2.618033988749895)
- middle:  u = v = 1.213411662762230 (real root of u^3 + u - 3 = 0)

Facts established here (keep consistent):

- With n = 2, subtracting the equations gives (u - v)(1 - uv) = 0; the
  asymmetric states have u + v = alpha, uv = 1, so u, v = (alpha ± sqrt(alpha^2 - 4))/2.
  Bistable (three steady states) iff alpha > 2; at alpha = 2 all merge at (1, 1)
  where det J = 0.
- Jacobian at the middle state: off-diagonals -1.191059, det J = -0.4186,
  eigenvalues 0.191 and -2.191 (a saddle; unstable).
- Jacobian at "U on": eigenvalues -1.6667 and -0.3333 (stable node; same at "V on").
- ODE basins: the diagonal u = v separates the two stable states (symmetry).
  Newton basins on [0, 4]^2: 21 % "U on", 57 % middle, 21 % "V on" (fig06-basins.svg).
- Newton from (3, 0.5): x1 = (2.6699267, 0.3594132); converges in 5 steps.
- Nearly singular start (1.9, 0.5): det J = -0.0299, first step to (-5.136, 4.425), ends in the middle state.
- With n = 2.5, the "U on" state is (2.955096, 0.187363).
