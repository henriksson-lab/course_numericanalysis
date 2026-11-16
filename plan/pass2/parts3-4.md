# Pass 2 report: parts3-4 (lu, pivoting-conditioning, structured-sparse, least-squares-qr, svd)

Every cross-reference in `plan/xref/<slug>.md` for the five units was checked against its target unit and tagged `[ok…]` / `[fixed: …]`. In `lec.qmd`, `exercise.qmd` and `quiz.qmd`, the first italic title mention per slide/page and in "Where this leads" is now a link `[*Title*](../<slug>/lec.qmd)`. No slides were added, no code was changed (except the display of one Rust tab) and no figures changed. `cargo build --examples` gives no warnings in these files.

## lu
- **Explore E-tasks no longer give away by-hand answers.** The intro paragraph no longer quotes the 13/15 operation counts, and it now tells students to do B3(b), B6(a) and B5 on paper first. E1 and E2 now ask students to *check* their own B3(b) counts against the demo. E5 now says to do B5(a)–(c) first, then check L, U and the diet with the swap button, and find the first step where they differ. The E5 answer adds the usual slip (forgetting to swap b).
- **C7 Rust tab** now shows the whole `rust/examples/lu_c7.rs`, helper functions included.
- Lecture notes ("Never solve(A) %*% b (2)") said κ(A4) "and why it matters are the next unit". pivoting-conditioning never uses A4. It now says what a condition number means is the next unit, which uses flow cytometry.
- Verified: the 1e9 ops/s rule and the "When timing misleads" slide (convergence-cost), NaN/Inf (floating-point), backward stability (error-conditioning). Also checked: stiffness reuses LU / modified Newton, and diffusion-implicit has the same tridiagonal matrix every step. B6(c) still treats that matrix as dense, as stated.

## pivoting-conditioning
- **The two-dye teaser matches error-conditioning exactly.** Same A = [1 0.95; 0.95 1], b = (1.95, 1.95) → (1.96, 1.94), x = (1, 1) → (1.2, 0.8), κ = 39 and the same R code; it also matches error-conditioning B9.
- The "Part 1" mention on the "condition number of a matrix" slide is now a link to error-conditioning.
- "Where this leads" for least-squares-qr now notes that unit teaches least squares on calibration curves, not spectra.
- Self-contained tabs: C2, C5 and C6 (R/Python) now point to `rel()` (C1) and `panel()` (C5).

## structured-sparse
- **Iterative-solvers preview verified.** iterative-solvers' recap quotes 4160 Jacobi iterations and 240 CG iterations on the 1000-cell morphogen, Jacobi diverging on [1 2; 2 1], and C8's `jacobi()`/`cg()`. All match. Gauss–Seidel, GMRES and preconditioning are covered there, as claimed.
- **Part 10 Thomas reuse verified.** diffusion-implicit (lecture and exercise) and advection-diffusion use `thomas(a, b, c, d)` with the same convention (a[1] and c[n] unused); the Rust signature is identical. steady-state-2d quotes the 30×30 fill-in numbers (4380 / 27 029 / 10 231) correctly.
- Fixed: the notes on the morphogen-matrix slide said diffusion-explicit gets "exactly the same matrix" and does a solve per time step. Now: the same tridiagonal pattern with the opposite sign and no degradation term, and diffusion-implicit does the per-step solve (1 + 2r on the diagonal).
- The C3 solution's "Part 10 uses it" now links to diffusion-implicit.
- Self-contained tabs: C4 and C5 (R/Python) now point to `thomas()` (C3).

## least-squares-qr
- The lecture notes and "Where this leads" implied nonlinear-least-squares fits the 4PL calibration curve. It fits a PK model, so the wording now says models like the 4PL are fitted by the methods of that unit.
- "Where this leads" bold-italic titles are now links. The Cholesky, κ recap, Runge, Gauss–Newton, R-based standard error and catastrophic-cancellation claims were all checked against their units and are true.

## svd
- All references were verified: least-squares-qr's `qr.solve`, the Rust `lstsq_qr` (same `q_tr_mul` + `solve_upper_triangular` pattern as `least_squares_qr_background.rs`), KKT in constrained-lp (which in turn points back to NNLS), the Jacobian's singular values in identifiability, irlba/Lanczos in eigen-algorithms, and the sister-course PCA remark (eigen-algorithms agrees). Only links were added.

## Problems needing changes outside this group
- None found. (diffusion-implicit calls Thomas every step instead of factoring once; it explains why, and that is consistent with lu B6(c)'s "treat as dense".)
