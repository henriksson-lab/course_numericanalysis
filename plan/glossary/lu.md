**Gaussian elimination** — Solving a linear system by subtracting multiples of one row from the rows below it until the matrix is upper triangular.
**Pivot** — The diagonal entry $u_{kk}$ used at step $k$ of elimination to zero out the entries below it.
**Multiplier** — The factor $\ell_{ik} = u_{ik}/u_{kk}$ by which the pivot row is multiplied before subtracting it from row $i$.
**LU factorisation** — Writing a square matrix as $A = LU$, with $L$ unit lower triangular (the multipliers) and $U$ upper triangular (the eliminated matrix); with row swaps, $A = PLU$.
**Lower / upper triangular matrix** — A matrix whose entries above (lower) or below (upper) the diagonal are all zero; unit lower triangular has 1s on the diagonal.
**Forward substitution** — Solving a lower triangular system $L\mathbf y = \mathbf b$ from the first row down, about $n^2$ operations.
**Back substitution** — Solving an upper triangular system $U\mathbf x = \mathbf y$ from the last row up, about $n^2$ operations.
**Factor once, solve many** — Computing the LU factorisation ($\approx \tfrac23 n^3$ operations) once and reusing it for every new right-hand side ($\approx 2n^2$ each).
**Permutation matrix** — An identity matrix with its rows reordered; multiplying by it swaps rows, and $P$ in $A = PLU$ records the row swaps of pivoting.
**Zero pivot** — A pivot equal to 0, which stops elimination without row swaps; it does not by itself mean the matrix is singular.
**Partial pivoting** — Swapping rows so that the pivot is the entry of largest absolute value in its column (treated in Pivoting, norms and conditioning).
**Residual** — The vector $\mathbf r = \mathbf b - A\mathbf x$ for a computed $\mathbf x$; a stable solver makes it as small as rounding allows.
**Isotope mixing model** — A linear system expressing a consumer's isotope signatures as a weighted average of its diet sources, with the weights (diet fractions) summing to 1.

## Install

No packages beyond those on `intro/install.qmd` (R: Matrix, shipped with R; Python: numpy, scipy; Rust: nalgebra, rand, rand_distr).
