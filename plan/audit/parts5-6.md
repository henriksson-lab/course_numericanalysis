# Pass 3 audit: parts5-6 (units 13–17)

Method: every R/Python block of lec/exercise/quiz was run cumulatively per page and language (`Rscript --vanilla`, `python3 -I`; R expression by expression so that intentional errors do not stop the page). Shown output comments were compared with the real output automatically, then by eye. Lecture snippets that are not self-contained were re-run with the setup rebuilt from the figure script or the exercise. Machine note: the load average was 200–270 on 40 cores for the whole run, so timings were inflated 5–10×.

## power-iteration (13)

- Code: lec, exercise (R and Python, Background + C1–C7) and quiz Q5 all match. Lec "shift exactly 1": the shown error line was simplified (`Error: system...`) → now `Error in solve.default(...): system is computationally singular: ...` (message text was already correct).
- Figures: `make-power-iteration-figs.R` runs clean; all 6 figures exist. fig06 note checked: the Rayleigh quotient reaches 1e-16 at step 24–25, the sum ratio ~48 ("about 50").
- By hand: re-did the lec 2×2 tables (1.75/1.4286/1.525/1.4918, eigen-decomposition 3/8, 1/8), inverse-iteration shift table, and B1–B6. All correct, with one fix: B2(c) said the error 7e-4 is "about 0.6·0.5^10" (= 5.9e-4) → "exactly 3.6·0.5^10/(5 + 0.5^10) ≈ 0.72·0.5^10".
- Quiz Q1–Q12: correct (0.9 rate → 22 steps/digit; Q9 rates; Q11 1.6e11 vs 5.6e9).
- Text vs numbers: rate 0.654 / 5.4 steps per digit / 55 steps, damping 1.53, 3.1 days, 20.6 days, λ2(Pslow) = 0.998, factor 500, inverse-iteration error table 3.2e-4 / 1.3e-7 / 5.1e-11, RQI table (re-run: 1.600000000000001 / 1.607466805519397 / 1.607467358251512, residuals 0.082 / 6.9e-4 / 4.4e-10), C eigenvalues 1.607/0.811/0.581, eigenvectors of P^T for 0.8 and 0.75, 8 TB, 5e9 vs 8e6. Fix: overflow step count in notes 4346 → 4347 (709.78/ln 1.177359).
- Audit notes: RQI overlap was handled in pass 2; Gupta et al. is a literature item; the 17 s vs 0.3 s timing is already labelled "one run ... your numbers will differ" in the notes.
- Rendered lec/exercise/quiz (retries after glob errors). No Rust changes.

## eigen-algorithms (14)

- **Recap slide (rewritten in pass 2)**: checked against power-iteration. λ1 = 1.177, rate 0.654, 55 steps, the semelparous fish cycling, the slow chain at 5406 steps, cubic RQI for symmetric A all match that unit. `power_rq(A, c(1,1,1))` re-run with power-iteration's own function gives 4.732051 in **51** steps. Rate 3/4.73 = 0.634.
- Code: lec re-run with the simulated X (exercise C6), `power_rq` and `ritz`. All outputs match (QR iterates, A^5 check, P stall, complex block, transit chain 0.1/0.1, kappa 1184131481, perturbation 0.09994035/0.09995391, Ritz values, eigs_sym/irlba 335.2194 ..., mprod 84/336), except the irlba residual 9.557529e-15: this reproduces with `OPENBLAS_NUM_THREADS=1` (as the notes say), and gives 7.7e-15/9.4e-15 multithreaded. Added a note that the last digits of this rounding-level number depend on the BLAS and its thread count.
- Exercise R/Python: all non-timing outputs match. The 1e-15 relative differences in C6 vary with BLAS threading (e.g. 1.53e-15 vs 1.36e-15), and the timings (C6, C7) vary with load. C6 text now says that rounding-level digits depend on the BLAS/threads as well; C7 text now says the printed times are from one machine.
- Figures: `make-eigen-algorithms-figs.R` runs clean; all 7 figures exist (fig07 uses typed-in timings, as documented).
- By hand B1–B6 and quiz Q1–Q12: all correct (B1 RQ = [[2.8, 0.6], [0.6, 1.2]], 25 more steps / 26 total = C2's 26; B4 ε^{1/n}; B5 Ritz 4.618/2.382 and 3.618/1.382; B6 memory/flops; Q2 83/13 steps; Q8 3.4e-5).
- Text vs numbers: 62/7/4 steps (C2 output), 7e-3 → 1e-7 → 6e-22, 0.025 rounding radius, 3.2 GB and half an hour, 1.6 GB, 2e7 flops/step, predator–prey period 8.908 and 0.64 per cycle. Fix: "λ2/λ1 = 0.99 needs about 3500 steps for 16 digits" → "about 3700" (ln 1e-16 / ln 0.99 = 3666).
- Audit notes: irlba/RSpectra are on intro/install.qmd; the timings typed into the fig script are documented in the script and in the notes; the recap overflow concern is cleared by the coordinator's overflow check (no overflow, slide unchanged in size).
- Rendered lec/exercise/quiz (retries after glob errors). No Rust changes.

## polynomial-interpolation (15)

- Code: all exercise outputs (C1–C7, R and Python) match. The n = 40 / n ≥ 30 entries differ between R and Python, which the page already describes as language-dependent rounding noise. The lec snippets were re-run with the fig-script setup and all match (od values, p(3), p(5), approx, Chebyshev, noisy fit 0.3937819/0.03161686, κ = 8.4e20, the `solve` refusal).
- Figures: the script runs clean (a harmless "y ≤ 0 dropped from log plot" warning); all 8 figures exist; pass 2 changed no figure data.
- By hand (Lagrange 0.4181, divided differences 0.4040, error bounds 0.037/0.086), B1–B6 and quiz Q1–Q11: all correct.
- Text vs numbers: all lec/exercise numbers match. Fixed 3 Explore answers in exercise.qmd: E1 range "0.03 and 0.5" → "0.03 and 0.4 (from n = 3 on)"; E2 "factor 1.6–1.8" → "about 1.6"; E4 Lebesgue constant "about 40" → "about 48".
- Audit note (figure slides not screenshotted) is covered by the coordinator's overflow check. Deck 23 slides. All 3 pages render. No Rust changes.

## splines (16)

- Code: running each page top to bottom exposed a session clash. The C1 and C3 solutions overwrote `x`/`y`, so C2 and C5 (which use the Background `t, y`) failed ("lengths differ" in R, broadcast errors in Python). Renamed C1's variables to `xh, yh` and C3's to `xk, yk` (R and Python). After that every shown output matches. The lec snippets were re-run with the fig-script/exercise setup and all match.
- Figures: the script runs clean; all 10 figures exist; the quoted figure numbers (fig03 slopes, fig06 df, fig07 9.57 vs 9.42, fig09, fig10) match.
- By hand B1–B7 and quiz Q1–Q12: correct.
- Text vs numbers, fixed: lec "A biology problem first" comment "true rate between 0.04 and 0.24" → "0.007 and 0.24" (N'(10) = 0.0073); lec "Boundary conditions" "about 5× larger" → "about 4×" (1.31e-3 vs 3.5e-4); lec notes GCV lag-1 correlation −0.18 → −0.19; exercise E1 shrink factors → 4.7, 3.7, 3.7, 3.0.
- Audit note (tradeSeq/Monocle) is a literature check, skipped. Deck 26 slides. All 3 pages render. No Rust changes.

## fourier-fft (17)

- Code: all non-timing outputs in lec (setup rebuilt from the fig script, seed 3), exercise (R and Python, C6 seeds match) and quiz Q9 match.
- Figures: the script runs clean and all 8 figures exist. fig01–03 and fig05–08 came out byte-identical. fig04 (timings) was regenerated with timings distorted by the load, so the previous copy was restored.
- By hand (n = 4 DFT, aliasing, Nyquist table), B1–B6, C4/C5 numbers and quiz Q1–Q12: correct.
- Text vs numbers, fixed: the trace ends near 26, not 20 (true signal at 119.5 h is 25.6). Fixed in the lec filtering-at-the-edges slide and in the exercise C6 text.
- Timing remarks added where they were missing: lec bad-`fft`-length slide (code comment and notes), exercise C2, quiz Q10 ("on our laptop").
- Audit note "R fft prime-length slowness": re-measured, 65537 vs 65536 samples is about 1600–10 000× under load, so "more than 1000×" holds. The period and sampling items are literature checks.
- Deck 25 slides; all 3 pages render. No Rust changes.

## Unresolved

- fourier-fft fig04 and the timing-slide claims (~1000×, about 1 s for n = 4×10^6, crossover near n = 256) could not be re-measured reliably under this load. Re-run `make-fourier-fft-figs.R` on an idle machine to settle them.
- Demo-derived numbers in the fourier-fft notes and its E1–E4 answers were not re-checked. The power-iteration and eigen-algorithms demo answers (E1–E5) were checked only where they coincide with code output.
