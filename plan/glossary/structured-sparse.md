**Symmetric positive definite (SPD)** — a symmetric matrix with x^T A x > 0 for every non-zero x; equivalently, all eigenvalues positive (e.g. covariance, kinship).
**Positive semidefinite** — symmetric with x^T A x >= 0; some eigenvalues may be zero (e.g. a covariance from fewer samples than variables).
**Cholesky factorisation** — A = L L^T with L lower triangular, for SPD A; half the cost of LU and no pivoting needed.
**Log-determinant** — log det A; from Cholesky as 2 * sum(log(diag(L))), avoiding the under- or overflow of det A itself.
**Shrinkage (ridge, nugget)** — adding lambda I to a covariance matrix to make it positive definite.
**Banded matrix** — a matrix whose non-zeros lie within p diagonals of the main diagonal; elimination without pivoting stays inside the band.
**Tridiagonal matrix** — a banded matrix with only the main diagonal and the diagonals directly above and below it.
**Thomas algorithm** — Gaussian elimination for a tridiagonal system without pivoting, about 8n operations.
**Diagonal dominance** — |a_ii| >= sum of |a_ij| over j != i in every row; makes elimination without pivoting and Jacobi iteration safe.
**Sparse matrix** — a matrix with so many zeros that storing only the non-zeros pays off.
**nnz** — the number of stored non-zero entries of a sparse matrix.
**Compressed sparse column (CSC)** — sparse storage by columns: values, row indices and one pointer per column (plus one); about 12 bytes per non-zero.
**Fill-in** — zero entries of A that become non-zero in its factors during elimination.
**Fill-reducing ordering** — a renumbering of the unknowns (e.g. AMD, COLAMD, nested dissection) chosen before factorising to keep fill-in small.
**Arrow matrix** — a sparse matrix with a full first (or last) row and column plus the diagonal, as from a hub gene linked to all others.
**Iterative solver** — a method that improves an approximate solution using only matrix-vector products, instead of factorising.
**Jacobi iteration** — x <- x + D^{-1}(b - A x) with D the diagonal of A; converges for diagonally dominant A.
**Conjugate gradient (CG)** — an iterative solver for SPD systems that chooses a new search direction each step; converges much faster than Jacobi.
**Residual** — r = b - A x; a cheap check of a computed solution and the stopping criterion of iterative solvers.

## Install

None. Matrix ships with R (already listed); scipy.sparse is part of scipy; sprs and nalgebra are in rust/Cargo.toml.
