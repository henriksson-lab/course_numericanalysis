**QR algorithm** — Repeat $A_k = Q_kR_k$, $A_{k+1} = R_kQ_k$; a sequence of orthogonal similarity transforms that converges to (block) triangular form with the eigenvalues on the diagonal.
**Similarity transform** — $A \mapsto S^{-1}AS$; it keeps the eigenvalues. With orthogonal $S$ it is $Q^TAQ$.
**Simultaneous (subspace) iteration** — Power iteration on several vectors at once, re-orthonormalised by QR after each step; the QR algorithm is this in disguise.
**Hessenberg matrix** — A matrix with zeros below the first subdiagonal; reducing to it first makes each QR step cost $O(n^2)$ instead of $O(n^3)$.
**Tridiagonal matrix** — Non-zeros only on the diagonal and the two neighbouring diagonals; the Hessenberg form of a symmetric matrix.
**Shift** — Running a QR step on $A - \mu I$ and adding $\mu I$ back, to speed up convergence near $\mu$.
**Wilkinson shift** — The eigenvalue of the trailing 2 × 2 block closest to $a_{nn}$; the shift used by LAPACK for symmetric matrices.
**Deflation** — Removing a converged eigenvalue (last row and column) and continuing on the smaller matrix.
**Real Schur form** — Upper triangular with 1 × 1 and 2 × 2 diagonal blocks; where the real QR algorithm stops on a nonsymmetric matrix, one 2 × 2 block per complex pair.
**Complex conjugate pair** — Eigenvalues $a \pm bi$ of a real matrix; in a linear ODE model an oscillation with period $2\pi/b$.
**Defective matrix** — A matrix with fewer independent eigenvectors than its size (e.g. a chain of identical compartments); its eigenvalues are extremely sensitive.
**Eigenvalue sensitivity** — How far eigenvalues move when the matrix entries change; at most $\|E\|$ for symmetric matrices, up to $\kappa(V)\|E\|$ (Bauer–Fike) or $\|E\|^{1/n}$ otherwise.
**Bauer–Fike bound** — For $A = VDV^{-1}$, a perturbation $E$ moves each eigenvalue by at most $\kappa(V)\|E\|$.
**Transit compartments** — A chain of identical first-order compartments used to model a delay (drug absorption, cell maturation).
**Krylov subspace** — $\text{span}\{\mathbf q, A\mathbf q, \dots, A^{m-1}\mathbf q\}$: all the vectors power iteration produces.
**Lanczos method** — Builds an orthonormal Krylov basis for a symmetric matrix with a three-term recurrence, giving a small tridiagonal matrix $T_m$.
**Ritz value** — An eigenvalue of $T_m = Q_m^TAQ_m$; an approximation to an eigenvalue of $A$ from the Krylov subspace.
**Ghost eigenvalue** — A spurious copy of a converged eigenvalue produced by Lanczos when its vectors lose orthogonality.
**Reorthogonalisation** — Removing from each new Lanczos vector its components along all earlier ones, to keep the basis orthonormal.
**Implicit restart** — Shrinking the Krylov basis while keeping its useful part, so that memory stays small (ARPACK, irlba).
**Matrix-free method** — A method that needs only a function computing $A\mathbf v$, never the matrix itself.
**Implicit centring** — Computing $X_c\mathbf v = X\mathbf v - \mathbf 1(\boldsymbol\mu^T\mathbf v)$ so that a sparse $X$ never has to be centred (and made dense).
**Truncated eigensolver** — A method that computes only the top $k$ eigenpairs or singular triplets (`irlba`, `RSpectra::eigs_sym`, `scipy.sparse.linalg.eigsh`, `svds`).
**Residual (eigenpair)** — $\|A\mathbf v - \lambda\mathbf v\|$; small means $(\lambda, \mathbf v)$ is an exact eigenpair of a nearby matrix.

## Install

- R: `irlba`, `RSpectra` (CRAN; not yet on `intro/install.qmd`). `Matrix` ships with R.
- Python: none beyond numpy/scipy (`scipy.sparse.linalg.eigsh`, `svds`, `LinearOperator`).
