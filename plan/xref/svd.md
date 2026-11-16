# Cross-references from `svd` (unit 12)

- svd/lec.qmd:24 (About this lecture) -> least-squares-qr (least squares by QR, assumes well-separated columns) [ok]
- svd/lec.qmd:74 (Least squares gives a negative cell type, code comment "previous unit") -> least-squares-qr (qr.solve for rectangular least squares) [ok]
- svd/lec.qmd:89 (notes, "as in the previous unit") -> least-squares-qr (qr.solve returns the least-squares solution) [ok]
- svd/lec.qmd:94 (Recap: conditioning of least squares) -> pivoting-conditioning (condition number, lose log10 kappa digits) [ok]
- svd/lec.qmd:94 (Recap: conditioning of least squares) -> least-squares-qr (normal equations square kappa; QR does not) [ok]
- svd/lec.qmd:177 (Why not take eigenvalues of A^T A?, "Same lesson as the normal equations") -> least-squares-qr (normal equations vs QR) [ok]
- svd/lec.qmd:400 (Proportions cannot be negative: NNLS, "QR, as in the previous unit") -> least-squares-qr (small least-squares solves by QR) [ok]
- svd/lec.qmd:404 (notes, KKT condition) -> constrained-lp (KKT / optimality conditions under constraints) [ok]
- svd/lec.qmd:450 (notes, Checklist) -> identifiability (same ideas for nonlinear models, Jacobian in the role of S) [ok]
- svd/lec.qmd:458 (Aside: PCA is an SVD) -> eigen-algorithms (irlba / Lanczos for top singular triplets of large sparse matrices) [ok]
- svd/lec.qmd:475 (Summary, Next) -> power-iteration (start of Part 5) [ok]
- svd/lec.qmd:460 (Aside: PCA is an SVD, "sister course") -> sister course quantbio-3MB00A `pca` (statistics of PCA; external, not a unit here) [ok]
- svd/exercise.qmd:300 (Hint 7) -> least-squares-qr (Rust NNLS built on the QR least-squares solve "from the previous unit") [ok]
- svd/exercise.qmd:555,696,849,1012,1148,1326 (Rust solutions C2–C7, comment in lstsq_qr "as in the previous unit") -> least-squares-qr (Householder QR least squares, q_tr_mul + solve_upper_triangular, same style as rust/examples/least_squares_qr_background.rs) [ok]
- svd/exercise.qmd:1442 (Where this leads) -> power-iteration (eigenvalue algorithms; singular values as sqrt of eigenvalues of A^T A) [ok]
- svd/exercise.qmd:1442 (Where this leads) -> eigen-algorithms (Lanczos / irlba for PCA of single-cell data) [ok]
- svd/exercise.qmd:1443 (Where this leads) -> identifiability (small singular values of the Jacobian = non-identifiable combinations) [ok]
- svd/exercise.qmd:1444 (Where this leads) -> constrained-lp (KKT conditions; NNLS optimality check is a special case) [ok]

## Running examples defined here

- **Bulk RNA-seq deconvolution signature** (used in lec, exercise, quiz, figures). 40 genes × 5 cell types
  (CD4T, CD8T, B, NK, Mono), built without random numbers so R, Python and Rust agree:
  start at 10; genes 1–16 = 100 for CD4T and CD8T; genes 17–24, 25–32, 33–40 = 100 for B, NK, Mono;
  multiply entry (g, j) by `1 + 0.5 * sin(g * j)` (g = 1..40, j = 1..5); set CD8T = CD4T, then
  CD8T gene 1 × 0.98 and gene 2 × 1.02. True proportions p = (0.30, 0.15, 0.15, 0.10, 0.30);
  bulk b0 = S p; "noisy" b = b0 * (1 + 0.05 * sin(37 g)). Key numbers: singular values
  642.93, 340.64, 278.84, 269.08, 2.8746; kappa = 223.66; v5 ≈ (1, −1, 0, 0, 0)/√2;
  least squares gives (0.543, −0.092, 0.149, 0.099, 0.299); NNLS gives (0.451, 0, 0.149, 0.099, 0.299);
  TSVD k = 4 gives (0.226, 0.226, 0.149, 0.099, 0.299); ridge λ = 10 gives (0.250, 0.201, 0.149, 0.099, 0.298).
  Exactly rank-deficient variant S0: CD8T column = CD4T column. Code: svd/lec.qmd slide
  "The simulated signature and bulk sample", Rust `signature()` in rust/examples/svd_setup.rs.
  Could be reused in `identifiability` (non-identifiable combinations) or `constrained-lp` (sum-to-one constraint).
- **2×2 example A = [[3, 0], [4, 5]]** with σ = √45, √5, v1 = (1,1)/√2, u1 = (1,3)/√10 (rotate–stretch–rotate figure, exercise B1).
