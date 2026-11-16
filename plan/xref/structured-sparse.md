# Cross-references from structured-sparse

- structured-sparse/lec.qmd:32 -> lu (previous unit, general dense LU)
- structured-sparse/lec.qmd:32 -> pivoting-conditioning (previous unit)
- structured-sparse/lec.qmd:98 -> pivoting-conditioning (element growth / why pivoting; Cholesky needs none)
- structured-sparse/lec.qmd:145 -> floating-point (underflow; sum of logs instead of product, applied to log det)
- structured-sparse/lec.qmd:194 -> pivoting-conditioning (ill-conditioning loses log10(kappa) digits)
- structured-sparse/lec.qmd:211 -> diffusion-explicit (same tridiagonal matrix from the diffusion equation)
- structured-sparse/lec.qmd:252 -> diffusion-implicit (Thomas function reused once per implicit time step)
- structured-sparse/lec.qmd:291 -> convergence-cost (reading slopes off log-log timing plots)
- structured-sparse/lec.qmd:309 -> pivoting-conditioning (failure of elimination without pivoting)
- structured-sparse/lec.qmd:"Recap: sparse storage" slide (l. 314, 323) -> convergence-cost (8 bytes/entry, CSC ~12 bytes/non-zero, centring destroys sparsity)
- structured-sparse/lec.qmd:406 -> lu (never compute the inverse to solve)
- structured-sparse/lec.qmd:443 -> steady-state-2d (2D five-point Laplacian, fill-in in a band)
- structured-sparse/lec.qmd:467 -> iterative-solvers (full treatment of Jacobi, Gauss-Seidel, CG, GMRES, preconditioning)
- structured-sparse/lec.qmd:497 -> pivoting-conditioning (small residual is not small error)
- structured-sparse/lec.qmd:497 -> convergence-cost (error vs iteration on log axes, linear convergence)
- structured-sparse/lec.qmd:508 -> iterative-solvers ("Which solver?" table)
- structured-sparse/lec.qmd:524 -> least-squares-qr (Next: normal equations, X^T X is SPD, Cholesky)
- structured-sparse/lec.qmd:527 -> eigen-algorithms (sparse PCA of single-cell data)
- structured-sparse/lec.qmd:527 -> iterative-solvers
- structured-sparse/lec.qmd:527 -> diffusion-implicit, steady-state-2d (Part 10, plain text "Part 10")
- structured-sparse/exercise.qmd:140 -> pivoting-conditioning (B3(c): partial pivoting)
- structured-sparse/exercise.qmd:629 -> diffusion-implicit (C3 solution: "Part 10 uses it", plain text)
- structured-sparse/exercise.qmd:871 -> floating-point (C5 solution: 2 lost to rounding next to 1e17)
- structured-sparse/exercise.qmd:1281 -> pivoting-conditioning (C8 solution: condition number, error vs residual)
- structured-sparse/exercise.qmd:1286 -> least-squares-qr (Where this leads)
- structured-sparse/exercise.qmd:1287 -> eigen-algorithms (Where this leads)
- structured-sparse/exercise.qmd:1288 -> diffusion-implicit (Where this leads)
- structured-sparse/exercise.qmd:1288 -> steady-state-2d (Where this leads)
- structured-sparse/exercise.qmd:1289 -> iterative-solvers (Where this leads)
- structured-sparse/quiz.qmd:54 -> floating-point (Q4 hint: underflow of a product)
- structured-sparse/quiz.qmd:209 -> iterative-solvers (Q11 answer)

## Running examples defined here

- **Morphogen in a row of cells** (tridiagonal system): n cells, each exchanges
  with its two neighbours at rate D = 1 per min and degrades at k = 0.01 per min
  (lifetime 100 min); production s = 1 in cell 1 only; no flux out of either end.
  Equations: -D c[i-1] + (2D + k) c[i] - D c[i+1] = s[i], with D + k on the
  diagonal of the two end cells. Decay length sqrt(D/k) = 10 cells; checks:
  residual, mass balance k * sum(c) = 1, c[11]/c[1] = 0.3680 (n = 100).
  Diagonally dominant and SPD.
- **Thomas function** `thomas(a, b, c, d)` (R; same signature in Python and
  Rust): all four vectors of length n, a = sub-diagonal with a[1] unused,
  b = diagonal, c = super-diagonal with c[n] unused, d = right-hand side.
  Python version copies b and d. Rust: `fn thomas(a: &[f64], b: &[f64], c: &[f64], d: &[f64]) -> Vec<f64>`
  in rust/examples/structured_sparse_c3.rs. Part 10 (diffusion-implicit) should reuse it.
- **LD matrix** Sigma_ij = 0.9^|i-j| for n = 2000 SNP markers: det underflows
  to 0, log det = (n-1) log(1 - 0.81) = -3319.802.
- **3 x 3 gene covariance** Sigma = [4 2 -2; 2 10 2; -2 2 6], L = [2 0 0; 1 3 0; -1 1 2], det 144.
- **Arrow (hub gene) matrix**: hub diagonal n, row/column of ones, others diagonal 2;
  n = 1000: chol has 500500 non-zeros hub first, 1999 hub last.
