# Cross-references from iterative-solvers

- iterative-solvers/lec.qmd:24 -> structured-sparse (glimpse of Jacobi and CG at the end of that unit)
- iterative-solvers/lec.qmd:32 -> structured-sparse, steady-state-2d (notes: prerequisites; 2D Laplacian rebuilt here in one slide)
- iterative-solvers/lec.qmd:37 -> structured-sparse (recap slide: morphogen row, Jacobi ~4160 vs CG ~240 iterations, [1 2; 2 1] diverges, exercise C8 jacobi()/cg())
- iterative-solvers/lec.qmd:59 -> steady-state-2d (notes: the culture-chamber oxygen problem, Michaelis-Menten version treated there)
- iterative-solvers/lec.qmd:82 -> steady-state-2d (notes: Kronecker construction of the five-point Laplacian)
- iterative-solvers/lec.qmd:106 -> structured-sparse (notes: CSC ~12 bytes per non-zero)
- iterative-solvers/lec.qmd:125 -> bisection-fixedpoint (notes: splitting = fixed-point iteration, |g'| < 1 analogue)
- iterative-solvers/lec.qmd:172 -> power-iteration (notes: repeated multiplication by T, residual ratio -> spectral radius)
- iterative-solvers/lec.qmd:172 -> eigen-algorithms (notes: dense eigen())
- iterative-solvers/lec.qmd:172 -> diffusion-explicit (notes: sine modes / eigenvalues of the Laplacian, cos(pi h))
- iterative-solvers/lec.qmd:185 -> pivoting-conditioning (condition number kappa)
- iterative-solvers/lec.qmd:190 -> diffusion-explicit (notes: same eigenvalues set the explicit time-step limit)
- iterative-solvers/lec.qmd:213 -> pivoting-conditioning (notes: small residual is not small error, error <= kappa * residual)
- iterative-solvers/lec.qmd:227 -> gradient-methods (notes: gradient descent zigzag on the quadratic; CG as its fix)
- iterative-solvers/lec.qmd:227 -> eigen-algorithms (notes: Lanczos = CG's recurrence)
- iterative-solvers/lec.qmd:238 -> convergence-cost (notes: linear convergence, slope on log plot)
- iterative-solvers/lec.qmd:273 -> structured-sparse (notes: the Cholesky algorithm that IC(0) truncates)
- iterative-solvers/lec.qmd:302 -> advection-diffusion (upwind scheme, Peclet number; perfused chamber)
- iterative-solvers/lec.qmd:307 -> advection-diffusion (notes)
- iterative-solvers/lec.qmd:333 -> least-squares-qr (notes: GMRES's small least-squares problem solved by QR)
- iterative-solvers/lec.qmd:333 -> eigen-algorithms (notes: Arnoldi / ARPACK)
- iterative-solvers/lec.qmd:344 -> fem-outlook (notes: algebraic multigrid on unstructured meshes)
- iterative-solvers/lec.qmd:377 -> steady-state-2d, reaction-diffusion, fem-outlook (Summary notes)
- iterative-solvers/exercise.qmd:21 -> structured-sparse (reuse jacobi()/cg() from exercise C8)
- iterative-solvers/exercise.qmd:219 -> advection-diffusion (C5: upwind differences)
- iterative-solvers/exercise.qmd:283 -> pivoting-conditioning (B4 solution: relative error <= kappa * relative residual)
- iterative-solvers/exercise.qmd:421 -> power-iteration (C1 solution: residual ratio = spectral radius)
- iterative-solvers/exercise.qmd:998 -> steady-state-2d, reaction-diffusion (Where this leads)
- iterative-solvers/exercise.qmd:999 -> fem-outlook (Where this leads)
- iterative-solvers/exercise.qmd:1000 -> eigen-algorithms (Where this leads: Lanczos, Arnoldi)
- iterative-solvers/quiz.qmd:126 -> structured-sparse (Q10 hint: fill-in of the full Cholesky factor)

## Running examples defined here

steady-state-2d/ did not exist when this unit was written, so the chamber was
built here; if steady-state-2d defines different numbers, one of the two
should be aligned.

- **Oxygen in a cell-laden gel (culture chamber)**: 1 mm x 1 mm square
  cross-section, O2 supplied through all four walls (c = 1, scaled by ambient
  c0 = 0.2 mM), uniform zero-order consumption. Deficit u = 1 - c solves
  -Laplace(u) = Phi, u = 0 on the walls, Phi = q L^2 / (D c0) = 12
  (D = 2000 um^2/s, q ~ 5 uM/s ~ 1e8 cells/mL x 50 amol/s). Grid: m x m
  interior points, h = 1/(m+1), equation multiplied by h^2:
  A = lap2(m) (4 on the diagonal, -1 neighbours, natural ordering i + m(j-1)),
  b = Phi h^2 * 1. Checks: centre c = 0.1161138 (m = 63), 0.1159863 (m = 127);
  14 % of the gel below c = 0.25 ("hypoxic") at m = 127.
  R: `lap1(m)`, `lap2(m) = kronecker(I, T) + kronecker(T, I)`, `chamber_rhs(m, Phi = 12)`.
- **Iteration counts (relative residual 1e-6)**, m = 15/31/63/127/255/511:
  Jacobi 705/2825/11302/45193; Gauss-Seidel 354/1414/5652/22598;
  CG 24/50/100/203/409/827; IC(0)-CG 13/23/40/73/144/275; multigrid V(2,2) 5 everywhere.
- **Aggregate chamber**: same, with D = 0.1 inside a disk of radius 0.3 mm at the
  centre (D = 1 elsewhere), A = -div(D grad) with D on faces. m = 63: CG 197,
  Jacobi-PCG 103, IC-CG 49.
- **Perfused chamber**: flow in +x, Pe = 50 (u = 100 um/s, L = 1 mm,
  D = 2000 um^2/s), upwind; m = 63: GMRES 116, GMRES(20) 187, CG diverges
  (true relative residual 314 after 300 iterations); centre c = 0.8800808.
