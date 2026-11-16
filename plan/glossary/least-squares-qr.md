**Overdetermined system** — A linear system with more equations than unknowns ($m > n$); in general no exact solution exists.
**Design matrix** — The matrix $A$ (R: `model.matrix`) with one row per measurement and one column per model parameter.
**Residual** — The part of the data the model misses, $\mathbf r = \mathbf b - A\mathbf x$.
**Least squares** — Choosing $\mathbf x$ to minimise the sum of squared residuals $\|\mathbf b - A\mathbf x\|^2$.
**Standard curve (calibration curve)** — A fit of instrument readings against known amounts, used to read off the amount in unknown samples.
**Column space** — The set of all vectors $A\mathbf x$; least squares projects the data onto it.
**Orthogonal projection** — The closest point in a subspace; the difference to it is perpendicular to the subspace ($A^T\mathbf r = \mathbf 0$).
**Normal equations** — $A^TA\mathbf x = A^T\mathbf b$, the square system whose solution is the least-squares fit; they square the condition number.
**Vandermonde matrix** — A matrix whose columns are powers $1, x, x^2, \dots$ of the data points; the design matrix of a polynomial fit.
**Orthonormal columns** — Columns that are mutually perpendicular and of length 1, so that $Q^TQ = I$.
**Orthogonal matrix** — A square matrix with $Q^TQ = I$; it rotates or reflects and keeps all lengths.
**QR factorisation** — $A = QR$ with $Q$ having orthonormal columns and $R$ upper triangular; least squares becomes $R\mathbf x = Q^T\mathbf b$.
**Gram–Schmidt** — Building orthonormal columns one at a time by subtracting the components along the earlier ones.
**Modified Gram–Schmidt** — Gram–Schmidt that projects the partly cleaned vector instead of the original column; loses orthogonality like $\kappa\varepsilon$ rather than $\kappa^2\varepsilon$.
**Loss of orthogonality** — How far a computed $Q$ is from orthonormal, measured by $\|Q^TQ - I\|$.
**Householder reflection** — The orthogonal matrix $I - 2\mathbf v\mathbf v^T/\mathbf v^T\mathbf v$ that mirrors a vector onto a coordinate axis; the building block of library QR.
**Läuchli matrix** — A classic example ($\delta$ entries below a row of ones) where forming $A^TA$ rounds away all the information.
**Aliased coefficient** — A coefficient `lm` cannot estimate because its column is a linear combination of others; reported as `NA`.
**Numerical rank** — The number of columns that are independent above a tolerance (`qr()$rank`, `np.linalg.matrix_rank`).
**Minimum-norm solution** — Among all least-squares solutions of a rank-deficient problem, the one with the smallest $\|\mathbf x\|$; what `np.linalg.lstsq` returns.
**Orthogonal polynomials** — Polynomial basis columns that are mutually perpendicular on the data (R: `poly()`); a well-conditioned alternative to powers of $x$.
**Amplification efficiency** — For qPCR, $E = 10^{-1/\text{slope}} - 1$ from the standard curve; 100 % means doubling per cycle.

## Install

No packages beyond those in `intro/install.qmd` (base R `qr`, `lm`, `poly`; numpy, scipy; nalgebra).
