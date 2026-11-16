# Pass 3 audit: parts3-4 (lu, pivoting-conditioning, structured-sparse, least-squares-qr, svd)

Method: every R block (lec, exercise, quiz) and every Python block (exercise, quiz) was extracted in page order and run cumulatively per page and language (`Rscript --vanilla`, `python3 -I`, one top-level expression at a time, so that an intentional error on a slide does not stop the rest). Printed output was compared line by line with the output comments. Each `make-<slug>-figs.R` was run on a scratch copy of the unit folder, and the resulting SVGs were compared byte for byte with the committed ones. All by-hand solutions and quiz answers were recomputed. Every number the text quotes from code or figures was checked against a rerun (figure data recomputed where needed). The machine was very busy during the audit (load 150–270 on 40 cores), so the timings measured here are not meaningful. Timings were therefore checked only for a "machine-dependent" remark. No Rust was changed, so no cargo rebuild was needed.

## lu
- Code: all R/Python outputs match, except timings (lec "Factor once", "Never solve(A) %*% b (1)", exercise C4/C5), which are machine-dependent.
- Fixed: the two timing slides now say in a code comment that the times are from our server, one thread, and will differ. The C5 solution says the same ("Measured (our server, one thread; your times will differ)"). C4 already said so.
- Fixed: "Never solve(A) %*% b (2)" said the inverse's residual is "$10^4$ times larger" than the LU residual, which is printed as exactly 0. It now says "at rounding level (here exactly 0); the inverse's is about $10^3$ times rounding level" ($4.7\times10^{-12}$ against $\varepsilon\max|A|\max|x| \approx 2\times10^{-15}$). Also checked: κ(A4) = 7.9e4 ("about 8e4").
- Figures: the script runs. fig01 and fig04 are unchanged. fig02 and fig03 are live timings and differ only because of machine load, so the committed versions were kept; the slide numbers (6.4 s, 0.044 s, 145 solves, 4 min vs 4 s) agree with the committed data.
- By hand: B1–B7, E1–E5 and Q1–Q11 recomputed. All correct, including the 13/15 operation counts, 23 h vs 2 min, 280 vs 163, ratio 9.74, and 202 diets on the boundary.
- Progress notes: the C7 Rust tab was already fixed in pass 2. The lu demo answer leaks were handled in pass 2 (students are told to do the tasks on paper first). lu/exercise.html was re-rendered.

## pivoting-conditioning
- Code: all R/Python outputs match. (The intentional `backsolve` error on "A zero pivot, a tiny pivot" is shown correctly.)
- Fixed: notes on "Detector noise in unmixed cells" gave the row norm of $A^{-1}$ as 14.2. It is 14.15, so now 14.1, matching exercise C7.
- Fixed: "Rule of thumb" slide: "16 digits survive κ = 10^6" (ambiguous, reads as all 16) → "Rounding alone still leaves 10 of 16 digits at κ = 10^6".
- Figures: the script runs, and all 7 SVGs are byte-identical. The κ table (1.6 … 7.2e8), sd 0.015/0.142 and κ(A1) ≈ 2.6 were checked.
- By hand: B1–B7 (3-digit arithmetic, PA = LU, norms, κ = 39, worst direction, residual vs error, zero pivot) and Q1–Q11 are all correct. The C4 claim that `kappa()` is "up to 17 % too high" matches the shown output (3.5 vs 3).
- Progress note "Cytek complexity index": a literature claim, left for factcheck.md.

## structured-sparse
- Code: outputs match, except:
  - Timings (lec "Sparse solves in R", exercise C4). Fixed: the slide code now says "our server; yours differ". C4 already said so.
  - Rounding-level digits that depend on the BLAS thread count: the `eigen(C)` near-zero eigenvalues on the lecture slide (match a default multi-threaded run), and numpy `slogdet`/Cholesky log-det last digits in C2. These are not errors. Added a remark to the C2 solution that the last printed digits depend on the BLAS library and thread count.
  - The "glimpse of iterative solvers" slide code says "built as before", so it is not self-contained. Rebuilt separately: Jacobi needs 4159 iterations, error 8.4e-10, Jacobi rate 0.995, residual of the 100-cell Thomas solve 1.8e-15. All as stated.
- Figures: the script runs, and all 5 SVGs are byte-identical (fig03 timings are typed in from an earlier measurement). Text numbers checked: 6.3 s, 0.55 s, 4380/27 029/10 231 (recomputed), and Jacobi 4158 / CG 237.
- By hand: B1–B7 (Cholesky, Thomas with fractions, CSC, fill-in, Jacobi, cost) and Q1–Q11 are all correct.

## least-squares-qr
- Code: all R/Python outputs match. (The intentional `solve` errors and the `qr.solve` error are shown as on the page.)
- Text numbers checked against reruns: sums 517.56/1790.19, RSS 1.2232 and 0.6266, σ̂ = 0.18, efficiency 99.5 %, 15 600 copies, κ cubic 5.8e9 → 72, 72/506/4624, κ(poly) = 4.24 at any degree, κ = 2e7 → errors 6.8e-3 vs 2.1e-10, `solve` refuses from degree 11 and `qr.solve` from degree 13.
- Figures: the script runs, and all 7 SVGs are byte-identical.
- By hand: B1–B6 and Q1–Q11 are all correct (normal equations, Läuchli, Gram–Schmidt, Householder, Vandermonde bounds, dependent column).
- Progress notes: overflow is covered by the brief ("no slide overflows"). xref line drift is not a correctness issue. The demo drag interaction could not be tested headless and is still untested.

## svd
- Code: all R/Python outputs match (the intentional `qr.solve` error is shown). The scipy 1.15.2 `nnls` bug quoted in C6 still reproduces with the installed scipy 1.15.2 (returns (1.988, 0.007) with residual 0; correct is (1.995, 0) with residual 0.00707).
- Fixed: B3(c) amplification "50/0.35 = 141" → "50 %/0.352 % = 142" (exact 142.1), and 0.35 % → 0.352 %.
- Fixed: B4(a) "C-code gives (0.997, 1.007)" (no such C task) → "the exact SVD gives (0.997, 1.007)" (checked).
- Figures: the script runs, and all 7 SVGs are byte-identical. fig07 statistics recomputed: sd 0.86/0.004/0.065/0.21, sum sd 0.007, NNLS zeros 37.8 %/37.4 %, TSVD 0.225. TSVD table, ridge sums (0.957 at λ = 100, 0.738 at 300, 0.232 at 1000) and σ5 = 6.6e-14 for the identical signature were also checked.
- By hand: B1–B6 and Q1–Q11 recomputed; correct apart from the two fixes above.
- Progress notes: the tool names (EPIC, quanTIseq, MuSiC, CIBERSORT) are literature claims, left for factcheck.md. nnls is already on the install page.

## Checks
- All 15 pages rendered with `quarto render` (some needed a retry because of concurrent-render temp-dir errors).
- No slides were added. Slide counts are lu 22, pivoting-conditioning 26, structured-sparse 26, least-squares-qr 27, svd 27.
- No Rust files were touched.

## Not resolved
- Timing figures (lu fig02/fig03) were not re-measured: the machine was overloaded, and the committed figures agree with the text.
- The least-squares-qr demo drag interaction is untested (needs a browser).
