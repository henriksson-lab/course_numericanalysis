# Cross-references from `lu` (Gaussian elimination as LU)

- lu/lec.qmd:24 (About this lecture) -> newton-systems (one linear solve J dx = -f per Newton step) [ok, linked]
- lu/lec.qmd:31,35 (About this lecture + notes) -> pivoting-conditioning (pivoting deferred; zero pivot only as teaser) [ok, linked]
- lu/lec.qmd:105 (notes, "Many consumers, one matrix") -> stiffness, diffusion-implicit, newton-systems (same matrix every time step / frozen Jacobian) [ok, linked]
- lu/lec.qmd:249,252 ("What libraries do: Matrix::lu" + notes) -> pivoting-conditioning (why the library reorders rows, partial pivoting) [ok, linked]
- lu/lec.qmd:271 (notes, "Counting the cost") -> convergence-cost (1e9 operations per second rule) [ok, linked]
- lu/lec.qmd:282 (notes, "Measured: slopes") -> convergence-cost ("When timing misleads" slide; doubling ratios, threading) [ok, linked]
- lu/lec.qmd:304 (notes, "Factor once, solve many") -> stiffness, diffusion-implicit (Parts 8 and 10: right-hand sides arrive one at a time) [ok, linked]
- lu/lec.qmd:335 (notes, "Never solve(A) %*% b (1)") -> structured-sparse (inverse of a banded/sparse matrix is dense; LU stays banded) [ok, linked]
- lu/lec.qmd:359 (notes, "Never solve(A) %*% b (2)") -> error-conditioning (backward stability), pivoting-conditioning (kappa of A4) [fixed: notes now say what a condition number means is the next unit, which uses flow cytometry, not A4]
- lu/lec.qmd:370 (notes, "The residual check") -> pivoting-conditioning (conditioning; small residual vs error) [ok, linked]
- lu/lec.qmd:394 (notes, zero pivot) -> floating-point (0/0 = NaN, x/0 = Inf silently) [ok, linked]
- lu/lec.qmd:412,415 (tiny pivot slide + notes) -> pivoting-conditioning (partial pivoting), floating-point (epsilon) [ok, linked]
- lu/lec.qmd:426,429 (Checks slide + notes) -> pivoting-conditioning (condition number; spectral unmixing in flow cytometry) [ok, linked]
- lu/lec.qmd:441,444 (Summary) -> pivoting-conditioning (next unit) [ok, linked]
- lu/exercise.qmd:35 (Background fact 6) -> pivoting-conditioning (partial pivoting) [ok, linked]
- lu/exercise.qmd:129 (B6(c)) -> diffusion-implicit (implicit tracer diffusion, Part 10) [ok: diffusion-implicit has the same matrix every step (uses Thomas, does not factor once); B6 says treat as dense, solution notes the band]
- lu/exercise.qmd:242 (Solution B6(c)) -> structured-sparse (banded matrices are cheaper) [ok, linked]
- lu/exercise.qmd:484 (Solution C1) -> floating-point (-18.4, 11.6 not exact in binary) [ok, linked]
- lu/exercise.qmd:595 (Solution C2) -> pivoting-conditioning (why the largest pivot) [ok, linked]
- lu/exercise.qmd:717 (Solution C3) -> pivoting-conditioning (condition number explains noise amplification in a 4-source mixing model) [ok, linked]
- lu/exercise.qmd:981 (Solution C6) -> error-conditioning (backward stability), pivoting-conditioning (ill-conditioning) [ok, linked]
- lu/exercise.qmd:1110 (Solution C7) -> pivoting-conditioning (partial pivoting) [ok, linked]
- lu/exercise.qmd:1115-1118 (Where this leads) -> pivoting-conditioning, structured-sparse, stiffness, diffusion-implicit, newton-systems, least-squares-qr [fixed: links]
- lu/quiz.qmd:82 (Q6) -> pivoting-conditioning (why libraries reorder rows) [ok, linked]

## Running examples defined here

- **Isotope mixing model, 3 sources × (δ13C, δ15N)** (lecture): rows = (sum to 1, δ13C, δ15N); sources C3 plants (−26, 2), maize/C4 (−10, 6), fish (−18, 14) ‰ (trophic-enrichment corrected).
  `A = rbind(c(1,1,1), c(-26,-10,-18), c(2,6,14))`; no-pivot LU: L = [1 0 0; −26 1 0; 2 1/4 1], U = [1 1 1; 0 16 8; 0 0 10].
  Consumers (columns of B): (1, −19.6, 5.6) → f = (0.5, 0.3, 0.2); (1, −14.8, 6.8) → (0.2, 0.6, 0.2); (1, −17.2, 11.2) → (0.1, 0.2, 0.7); (1, −24, 9) → (0.55, −0.20, 0.65), outside the mixing triangle.
- **Zero-pivot variant**: C3 plants (−26, 2), freshwater fish (−26, 10), maize (−10, 6): `rbind(c(1,1,1), c(-26,-26,-10), c(2,10,6))`; zero pivot at step 2; consumer (1, −20, 6.5) → (0.25, 0.375, 0.375).
- **Near-duplicate 4-source / 3-isotope matrix** A4 (δ13C, δ15N, δ34S): columns (−26, 2, 4), (−26+d, 2+d, 4+2d), (−10, 6, 3), (−18, 14, 19), f = (0.3, 0.2, 0.4, 0.1); d = 0.01 in the lecture (kappa ≈ 7.9e4). Natural to reuse in pivoting-conditioning for the condition number.
- **Salt-marsh mixing model** (exercise): Spartina (δ13C −13, δ34S +6), phytoplankton (−21, +18), upland C3 (−27, +2); `rbind(c(1,1,1), c(-13,-21,-27), c(6,18,2))`; L = [1 0 0; −13 1 0; 6 −1.5 1], U = [1 1 1; 0 −8 −14; 0 0 −25]; mussels (1, −18.4, 11.6) → (0.4, 0.5, 0.1) and (1, −22.4, 7.6) → (0.2, 0.3, 0.5). 4-source version adds δ15N (6, 8, 4) and benthic microalgae (−17, 5, 10).

## Pass 2 (other fixes)

- [fixed: Explore intro no longer quotes the 13/15 counts; E1, E2 ask students to check their B3(b)/B6(a) counts; E5 asks them to do B5(a)–(c) first and check L, U and the diet with the swap button]
- [fixed: C7 Rust tab now shows the whole rust/examples/lu_c7.rs, helpers included]
