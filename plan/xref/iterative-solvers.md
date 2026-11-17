# Cross-references from iterative-solvers

Unit moved from Part 11 (optional) to Part 10 (core), between advection-diffusion
and reaction-diffusion. Recap and "About this lecture" now lean on
steady-state-2d (Kronecker Laplacian, sparse Cholesky fill-in table: nnz(A) 252 649,
nnz(L) 1 475 954 at h = 12.5 µm) and advection-diffusion (upwind weights 1 + P,
marching instead of a 2D solve); both verified against those lectures.
Next unit: reaction-diffusion; fem-outlook referred to as Part 11 (optional).

- iterative-solvers/lec.qmd:24 -> steady-state-2d, advection-diffusion (About this lecture {.smaller}) [linked]
- iterative-solvers/lec.qmd:33 -> structured-sparse, steady-state-2d, advection-diffusion (About this lecture {.smaller}) [linked]
- iterative-solvers/lec.qmd:38 -> structured-sparse (Recap: what we already have {.smaller}) [linked]
- iterative-solvers/lec.qmd:39 -> steady-state-2d (Recap: what we already have {.smaller}) [linked]
- iterative-solvers/lec.qmd:40 -> advection-diffusion (Recap: what we already have {.smaller}) [linked]
- iterative-solvers/lec.qmd:58 -> steady-state-2d (Biology: oxygen in a culture chamber {.smaller}) [linked]
- iterative-solvers/lec.qmd:81 -> steady-state-2d (The matrix: a 2D Laplacian {.smaller}) [linked]
- iterative-solvers/lec.qmd:124 -> bisection-fixedpoint (Jacobi and Gauss–Seidel as splittings {.smaller}) [linked]
- iterative-solvers/lec.qmd:171 -> power-iteration, eigen-algorithms, diffusion-explicit (The rate is a spectral radius {.smaller}) [linked]
- iterative-solvers/lec.qmd:184 -> pivoting-conditioning (Why they slow down as the grid is refined {.smaller}) [linked]
- iterative-solvers/lec.qmd:189 -> diffusion-explicit (Why they slow down as the grid is refined {.smaller}) [linked]
- iterative-solvers/lec.qmd:212 -> pivoting-conditioning (Failure: a small step is not a small error {.smaller}) [linked]
- iterative-solvers/lec.qmd:226 -> gradient-methods, eigen-algorithms (Conjugate gradient: the Krylov idea {.smaller}) [linked]
- iterative-solvers/lec.qmd:237 -> convergence-cost (CG against the splittings {.smaller}) [linked]
- iterative-solvers/lec.qmd:272 -> structured-sparse (Preconditioning {.smaller}) [linked]
- iterative-solvers/lec.qmd:301 -> advection-diffusion (Flow through the chamber: $A$ is not symmetric {.smaller}) [linked]
- iterative-solvers/lec.qmd:306 -> advection-diffusion (Flow through the chamber: $A$ is not symmetric {.smaller}) [linked]
- iterative-solvers/lec.qmd:332 -> least-squares-qr, eigen-algorithms (GMRES in one slide {.smaller}) [linked]
- iterative-solvers/lec.qmd:389 -> fem-outlook (Multigrid in one picture {.smaller}) [linked]
- iterative-solvers/lec.qmd:422 -> reaction-diffusion (Summary {.smaller}) [linked]
- iterative-solvers/lec.qmd:425 -> steady-state-2d, reaction-diffusion, fem-outlook (Summary {.smaller}) [linked]
- iterative-solvers/exercise.qmd:23 -> structured-sparse (Goal) [linked]
- iterative-solvers/exercise.qmd:27 -> steady-state-2d (Background) [linked]
- iterative-solvers/exercise.qmd:241 -> advection-diffusion (In code) [linked]
- iterative-solvers/exercise.qmd:332 -> pivoting-conditioning (Solutions) [linked]
- iterative-solvers/exercise.qmd:487 -> power-iteration (Solutions) [linked]
- iterative-solvers/exercise.qmd:1411 -> reaction-diffusion (Where this leads) [linked]
- iterative-solvers/exercise.qmd:1412 -> steady-state-2d (Where this leads) [linked]
- iterative-solvers/exercise.qmd:1413 -> fem-outlook (Where this leads) [linked]
- iterative-solvers/exercise.qmd:1414 -> eigen-algorithms (Where this leads) [linked]
- iterative-solvers/quiz.qmd:126 -> structured-sparse (None) [linked]

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
  Restarted: GMRES(10) 147, GMRES(20) 187, GMRES(50) 284; BiCGSTAB 88 steps
  (R; Python 86, Rust 88).
- **Signal chamber** (added when the unit became core, Part 10): the perfused
  chamber with a linearised positive feedback, A - sigma h^2 I, sigma = 100
  (same right-hand side). m = 63: GMRES(10), GMRES(20) stall at relative
  residual 0.976 / 0.973 after 3000 iterations; GMRES(50) 800; full GMRES 123;
  BiCGSTAB 103 (R; Python 100, Rust 101). Smooth bump v = sin(pi x) sin(pi y):
  v'Av/v'v = 0.005759 (perfused), -0.018655 (signal). sigma = 200: BiCGSTAB
  fails, full GMRES 133.

Pass 2: the gel chamber stays (different problem from steady-state-2d's lid chamber); the difference is now stated in lec notes (Biology slide) and in the exercise Background.
