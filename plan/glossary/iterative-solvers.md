**Iterative solver** — a method that improves a guess for the solution of Ax = b using only products with A, instead of factorising A.
**Splitting** — writing A = M - N with M easy to solve, giving the iteration x <- x + M^{-1}(b - Ax).
**Jacobi method** — the splitting with M = diagonal of A: every unknown is updated from the old values of its neighbours.
**Gauss–Seidel method** — the splitting with M = lower triangle of A: unknowns are updated in turn, using new values as soon as they exist.
**Iteration matrix** — T = I - M^{-1}A; the error of a splitting is multiplied by T in every step.
**Spectral radius** — the largest absolute value of the eigenvalues of a matrix, rho(T); a splitting converges if and only if rho(T) < 1, and gains one digit every ln 10 / (-ln rho) steps.
**Krylov space** — span{b, Ab, ..., A^{k-1} b}: all vectors that can be made from b with k - 1 products with A.
**Conjugate gradient (CG)** — Krylov method for symmetric positive definite matrices that picks the best point in the Krylov space with a short recurrence; needs about sqrt(kappa) iterations.
**A-conjugate** — two directions p, q with p^T A q = 0; CG's search directions are mutually A-conjugate.
**Preconditioner** — a matrix M close to A that is cheap to solve with; the solver works on M^{-1}A, which has a smaller condition number.
**Jacobi (diagonal) preconditioner** — M = diagonal of A; removes bad scaling between rows, nothing more.
**Incomplete Cholesky (IC(0))** — Cholesky factorisation in which every fill-in entry is dropped, so L keeps the sparsity pattern of A and LL^T only approximates A.
**GMRES** — generalised minimal residual method: the Krylov method for nonsymmetric matrices that minimises the residual; stores one basis vector per iteration.
**Restarted GMRES, GMRES(k)** — GMRES restarted from the current iterate every k steps to limit memory; can converge more slowly or stall.
**Stall (of restarted GMRES)** — the residual stops decreasing because each restart throws away the directions that were needed; possible when v^T A v <= 0 for some v.
**Symmetric part** — (A + A^T)/2; v^T A v = v^T (symmetric part) v. If it is positive definite, GMRES(k) converges for every k.
**BiCGSTAB** — bi-conjugate gradient stabilised: Krylov method for nonsymmetric matrices with short recurrences (fixed memory, two products with A per step) and an erratic residual; can stagnate or break down.
**Arnoldi process** — Gram–Schmidt orthogonalisation of the Krylov vectors; produces the basis and the small Hessenberg matrix used by GMRES.
**Smoother** — a few sweeps of Gauss–Seidel or damped Jacobi, which remove high-frequency (wiggly) error quickly but smooth error slowly.
**Multigrid** — solver that smooths on the fine grid and corrects the remaining smooth error on coarser grids (V-cycle); iteration count independent of the grid spacing.
**V-cycle** — one multigrid iteration: smooth and restrict down to the coarsest grid, solve there, interpolate and smooth back up.
**Algebraic multigrid (AMG)** — multigrid that builds its coarse levels from the matrix alone, without a geometric grid.

## Install

None for the exercise (R: Matrix; Python: numpy/scipy; Rust: sprs, nalgebra). Mentioned but not needed (and not installed): R packages `Rlinsolve`, `cPCG` (iterative solvers, incl. BiCGSTAB); Python `pyamg` (algebraic multigrid). `pracma` (installed; `gmres`, dense only). Python's `scipy.sparse.linalg.bicgstab` is used in C7.
