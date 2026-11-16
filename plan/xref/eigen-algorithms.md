# Cross-references from eigen-algorithms

- eigen-algorithms/lec.qmd:24 -> power-iteration (About: one eigenvalue at a time by repeated multiplication)
- eigen-algorithms/lec.qmd:74 -> sister course quantbio-3MB00A pca / pca-r (PCA statistics: number of PCs, scaling, interpretation)
- eigen-algorithms/lec.qmd:77 -> svd (notes: "Aside: PCA is an SVD" slide)
- eigen-algorithms/lec.qmd:99 -> power-iteration (notes: recap slide, power method, rate |lambda2/lambda1|)
- eigen-algorithms/lec.qmd:125 -> least-squares-qr (notes: Householder QR)
- eigen-algorithms/lec.qmd:184 -> power-iteration (notes: power iteration fails on a periodic Markov chain)
- eigen-algorithms/lec.qmd:207 -> power-iteration (notes: shifted inverse iteration, Rayleigh quotient)
- eigen-algorithms/lec.qmd:227 -> convergence-cost (notes: O(n^3) cost; "Part 1")
- eigen-algorithms/lec.qmd:262 -> phase-plane (notes: classification of Jacobian eigenvalues)
- eigen-algorithms/lec.qmd:322 -> pivoting-conditioning (notes: condition number kappa(V) in Bauer-Fike)
- eigen-algorithms/lec.qmd:357 -> iterative-solvers (notes: CG and GMRES use Krylov subspaces)
- eigen-algorithms/lec.qmd:416 -> svd (notes: irlba works on X, not X^T X)
- eigen-algorithms/lec.qmd:452 -> least-squares-qr, svd (notes: covariance route squares the condition number)
- eigen-algorithms/lec.qmd:474 -> error-conditioning (notes: backward error)
- eigen-algorithms/lec.qmd:499 -> polynomial-interpolation (Next: slide)
- eigen-algorithms/exercise.qmd:Where this leads -> phase-plane (Jacobian eigenvalues classify steady states)
- eigen-algorithms/exercise.qmd:Where this leads -> stiffness (spread of Jacobian eigenvalues)
- eigen-algorithms/exercise.qmd:Where this leads -> iterative-solvers (CG/GMRES, Lanczos/Arnoldi for Ax = b)
- eigen-algorithms/exercise.qmd:Where this leads -> svd (truncated SVD as best low-rank approximation); sister course (PCA statistics)

## Running examples defined here

- **3 x 3 symmetric test matrix** `A <- matrix(c(4,1,0, 1,3,1, 0,1,2), 3)`, eigenvalues 3 + sqrt 3 = 4.732051, 3, 3 - sqrt 3 = 1.267949 (already tridiagonal). Unshifted QR: 62 steps to deflate all; Rayleigh shift 7; Wilkinson 4.
- **Predator-prey Jacobian** `J <- rbind(c(-0.1, -1), c(0.5, 0))`, eigenvalues -0.05 +- 0.7053368i, period 8.908 (useful for `phase-plane`).
- **Transit-compartment chain** n = 10, rate 1: `K <- diag(-1, n); K[cbind(2:n, 1:(n-1))] <- 1`; feedback eps in K[1, n] moves all eigenvalues by eps^(1/n) (useful for `phase-plane`, `stiffness`, pharmacokinetics).
- **Simulated single-cell matrix** 20 000 cells x 2 000 genes, 10 types, `simulate_cells()` in `eigen-algorithms/make-eigen-algorithms-figs.R` / exercise C6 (seed 1): 13 % non-zero, top singular values of centred X 335.2194, 294.3642, ...; gap after 10 to ~72.4.
