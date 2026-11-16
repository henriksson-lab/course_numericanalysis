**Singular value decomposition (SVD)** — Factorisation A = UΣVᵀ of any matrix into orthonormal U and V and a diagonal Σ of non-negative singular values.
**Singular value** — A diagonal entry σᵢ of Σ; the stretch factor of A along the direction vᵢ (A vᵢ = σᵢ uᵢ).
**Left / right singular vectors** — The columns uᵢ of U (data space) and vᵢ of V (parameter space).
**Thin (economy) SVD** — The SVD with U of size m×n instead of m×m; what `svd()` and `np.linalg.svd(full_matrices=False)` return.
**Spectral norm (2-norm of a matrix)** — The largest stretch ‖A‖₂ = σ₁.
**Condition number (2-norm)** — κ(A) = σ₁/σₙ, the ratio of largest to smallest singular value.
**Numerical rank** — The number of singular values above a tolerance (max(m,n)·ε·σ₁ for rounding, the noise level for data).
**Rank deficiency** — Columns of a matrix that are (exactly or nearly) linear combinations of each other, so some parameter combinations are not determined.
**Collinearity** — Rank deficiency among predictors in a regression; shows as `NA` coefficients in `lm` or as small singular values.
**Pseudoinverse (Moore–Penrose)** — A⁺ = VΣ⁺Uᵀ, inverting only the non-zero singular values; A⁺b is the minimum-norm least-squares solution.
**Minimum-norm solution** — Among all least-squares solutions, the one with the smallest ‖x‖.
**Truncated SVD (TSVD)** — Least-squares solution keeping only the k largest singular values.
**Regularisation** — Changing an ill-conditioned problem to damp the noise-dominated components, at the price of bias.
**Tikhonov (ridge) regularisation** — Minimising ‖Ax − b‖² + λ²‖x‖²; damps components with σᵢ below λ.
**Filter factor** — The weight fᵢ multiplying SVD component i of a regularised solution; σᵢ²/(σᵢ² + λ²) for ridge, 1 or 0 for TSVD.
**Picard plot** — Log plot of σᵢ and |uᵢᵀb| against i; shows where noise takes over the data coefficients.
**L-curve** — Log–log plot of solution norm against residual norm over λ; its corner is a common choice of λ.
**Discrepancy principle** — Choose the regularisation so that the residual is about as large as the expected noise.
**Non-negative least squares (NNLS)** — Least squares subject to x ≥ 0.
**Active-set method (Lawson–Hanson)** — NNLS algorithm that holds some variables at 0 and solves ordinary least squares on the others, updating the sets until the optimality conditions hold.
**Signature matrix** — Genes × cell types matrix of reference expression profiles used to deconvolve bulk samples.
**Deconvolution (bulk RNA-seq)** — Estimating cell-type proportions p from a bulk profile b ≈ S p.

## Install

- R: `nnls` (non-negative least squares, `install.packages("nnls")`). `MASS` (for `MASS::ginv`) ships with R. `irlba` is mentioned on a slide but not used in code.
- Python: nothing beyond numpy, scipy, matplotlib (`scipy.optimize.nnls`).
