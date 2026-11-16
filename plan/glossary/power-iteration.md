**Eigenvalue / eigenvector** — A number $\lambda$ and nonzero vector $\mathbf v$ with $A\mathbf v = \lambda\mathbf v$: a direction the matrix only stretches.
**Dominant eigenvalue** — The eigenvalue of largest modulus, $|\lambda_1|$; strictly dominant if no other has the same modulus.
**Characteristic polynomial** — $\det(A - \lambda I)$; its roots are the eigenvalues (practical by hand only for 2 × 2 and 3 × 3).
**Leslie matrix** — An age-structured population matrix with fertilities in the first row and survival probabilities on the subdiagonal; $\mathbf n_{t+1} = L\mathbf n_t$.
**Long-term growth rate** — The dominant eigenvalue $\lambda_1$ of a population projection matrix; the population grows ($\lambda_1 > 1$) or declines ($\lambda_1 < 1$) by this factor per time step in the long run.
**Stable age (stage) distribution** — The dominant eigenvector of a projection matrix, rescaled to proportions; the age mix the population settles into.
**Semelparous / iteroparous** — Reproducing once and then dying / reproducing in several seasons; semelparity can make all eigenvalues of a Leslie matrix equal in modulus.
**Damping ratio** — $\lambda_1/|\lambda_2|$ of a projection matrix: how fast a population forgets its initial structure; its reciprocal is the convergence rate of power iteration.
**Transition matrix (stochastic matrix)** — A Markov chain matrix with $P_{ij}$ = probability of moving from state $i$ to state $j$; rows sum to 1.
**Stationary distribution** — Proportions $\boldsymbol\pi$ with $\boldsymbol\pi P = \boldsymbol\pi$: the eigenvector of $P^T$ for eigenvalue 1, rescaled to sum 1.
**Relaxation time** — How many steps a Markov chain needs to approach its stationary distribution, set by $|\lambda_2|$ (e.g. halving time $\log 0.5/\log|\lambda_2|$).
**Perron–Frobenius theorem** — A nonnegative matrix some power of which is entirely positive has a real, positive, strictly dominant eigenvalue with a positive eigenvector.
**Power iteration (power method)** — Repeating $\mathbf x \leftarrow A\mathbf x/\|A\mathbf x\|$; converges to the dominant eigenvector at rate $|\lambda_2/\lambda_1|$.
**Normalisation** — Rescaling the iterate every step (by its sum, 2-norm or largest entry) to avoid overflow or underflow; it does not change the direction.
**Eigen-residual** — $\|A\mathbf x - \lambda\mathbf x\|$; zero for an exact eigenpair, the quantity to stop on.
**False convergence** — A tiny change per step while the error is still large, typical when the convergence rate is close to 1; error $\approx$ change $\cdot\,c/(1-c)$.
**Rayleigh quotient** — $R(\mathbf x) = \mathbf x^TA\mathbf x/\mathbf x^T\mathbf x$, the best eigenvalue estimate for a given vector; error of order $\delta^2$ for symmetric $A$.
**Inverse iteration** — Power iteration with $(A - \sigma I)^{-1}$ (one linear solve per step); finds the eigenvalue nearest the shift $\sigma$.
**Shift** — The number $\sigma$ subtracted from the diagonal in inverse iteration; the closer it is to the wanted eigenvalue (and the further from the others), the faster the convergence.
**Rayleigh quotient iteration** — Inverse iteration with the shift updated to the current Rayleigh quotient each step; cubic convergence for symmetric matrices.
**Integral projection model (IPM)** — A population model with a continuous trait (e.g. body size) discretised into many classes; a large Leslie-like matrix.

## Install

No packages beyond those in `intro/install.qmd` (base R `eigen`, `solve`; `Matrix::lu` for factor-once inverse iteration; numpy, scipy.linalg `lu_factor`/`lu_solve`; nalgebra).
