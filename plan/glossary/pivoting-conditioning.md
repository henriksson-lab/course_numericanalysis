**Pivot** — the diagonal entry $a_{kk}$ that elimination divides by when clearing column $k$.
**Partial pivoting** — swapping rows before each elimination step so the largest entry (in absolute value) of the current column becomes the pivot; keeps all multipliers at most 1 in size.
**Permutation matrix** — the identity matrix with its rows reordered; multiplying by it reorders rows, as in $PA = LU$.
**Vector norm** — a measure of the length of a vector; the common ones are the 1-norm (sum of absolute values), 2-norm (Euclidean length) and $\infty$-norm (largest absolute value).
**Matrix norm** — the largest factor by which a matrix can stretch a vector, $\max \|A\mathbf v\|/\|\mathbf v\|$; the 1-norm is the largest column sum, the $\infty$-norm the largest row sum, the 2-norm the largest singular value.
**Singular value** — a stretch factor of a matrix along one of its principal directions; $\sigma_{\max}$ and $\sigma_{\min}$ are the largest and smallest.
**Condition number of a matrix** — $\kappa(A) = \|A\|\,\|A^{-1}\| = \sigma_{\max}/\sigma_{\min}$; the largest factor by which relative errors in $\mathbf b$ are amplified in the solution of $A\mathbf x = \mathbf b$.
**Ill-conditioned matrix** — a matrix with large $\kappa$; solving with it loses about $\log_{10}\kappa$ digits of the data's accuracy.
**Residual** — $\mathbf r = \mathbf b - A\hat{\mathbf x}$, how far a computed solution is from satisfying the equations; computable, unlike the error.
**Backward stable (linear solve)** — the computed solution is the exact solution of a system with $A$ changed by about machine epsilon; true of LU with partial pivoting in practice.
**Spectral unmixing** — recovering the amounts of several dyes from detector readings by solving $A\mathbf x = \mathbf b$, where $A$ holds the dyes' signals in each detector.
**Spillover spreading** — the extra spread of unmixed cytometry data caused by detector noise amplified through an ill-conditioned unmixing matrix.

## Install

None beyond `intro/install.qmd` (uses base R plus the `Matrix` package; numpy and scipy; Rust crates nalgebra, rand, rand_distr).
