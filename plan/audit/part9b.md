# Pass 3 audit: part9b (units 27–30)

Method: every R block (lec, exercise, quiz) and every Python block (exercise, quiz) was extracted and run cumulatively per page and language (`Rscript --vanilla`, `python3 -I`). Printed output was compared with the shown output in both directions. Lecture fragments that rely on objects defined off-slide were run with those objects supplied: `pstar`, `sse_min` and `sse_h` in gradient-methods, the `try()`-wrapped `nls` errors in nonlinear-least-squares and identifiability. The identifiability quiz was run after the exercise's setup block, as its text says. All four figure scripts run without errors. Every figure that lec/exercise reference exists, and all 28 SVGs come out byte-identical to the ones before the audit. All by-hand solutions and quiz answers were re-derived, and the numbers quoted in the prose were checked against the outputs. Rust was not touched.

## gradient-methods

- Code: all R/Python output matches, including the lecture with `pstar` (BFGS optimum) and the growth-curve `sse_min`/`sse_h` supplied.
- Fixed in lec: on the slide "Steepest descent with a fixed step", the comment said "NaN after 11 steps". In fact EC50 < 0, and so SSE = NaN, after 10 steps; the comment now says "EC50 < 0 after 10 steps".
- Fixed in lec: in the Newton table, SSE at steps 4–5 was 314.92395940. R and Python give 314.92395944.
- Fixed in exercise: C2 said "about 9 times more steps for every factor 10 in κ". The real counts are 74 → 840 → 9243, so it now says "about 11 times".
- Fixed in exercise: B2(c) said the Newton step goes "towards the saddle". It goes to x = −1.25, so it now says "towards (and past) the saddle".
- Checked, no change: B1–B8 arithmetic. Quiz Q4 (62 → ~6200 steps). The 3800-step estimate for κ = 550. 2/λmax = 5.4e-4. Pure Newton ends at EC50 86, hill 16, SSE 21188. Modified Newton takes 9 steps. In nM: κ = 1.9e6, 109 calls; with parscale 36 calls. Growth curve: κ 1.2e4 → 1.0e3. C1/C3/C5–C8 prose against the R and Python outputs (e.g. scipy NM fails from 7 starts, 6 of them reported as success).
- Progress-notes item: the "design demo typecheck error in demos/src/design/math.ts:109" (the design demo belongs to identifiability) is gone. `tsc --noEmit -p demos/tsconfig.json` passes cleanly.

## nonlinear-least-squares

- Code: all R/Python output matches. The lecture's unlabeled Gauss–Newton output fence matches too. `minpack.lm` is installed.
- Fixed in exercise: the C1 solution claimed that the 1.6e-9 difference from the ka = ke limit is "mostly the genuine effect". A 40-digit computation shows the true relative effect at t = 1 h is 4.5e-9. The rest, about −3e-9, is cancellation error. The sentence now says this.
- Fixed in exercise: C6(b) said the solver error at tol 1e-10 "is about 1e-10". R gives 5e-10 and scipy 3e-11, so it now says "below 1e-9 (R 5e-10, scipy 3e-11)".
- Added (optional pass-2 item, no new slide): one sentence in the exercise intro that the same methods and code fit the logistic growth curve (intro) and an SIR model (euler), with links. Also one sentence in the notes of the slide "Fitting an ODE model".
- Checked, no change: the by-hand GN table (8.11 → 0.285, steps 0.5021, 0.504413, 0.5044173). B1–B7 arithmetic. κ(J) 170, κ(JᵀJ) 2.9e4, log κ 8.7. The bad GN step (ka −2.22, V 44 mL, RSS 6.7e29, κ ≈ 2e12). LM: 15 passes, rejected 1, 2, 4, 13, 14. Log GN runs off to ka 2e-4, ke 617. Multistart 14/21/1 with RSS 64.5 = Σy². Solver-error figure: 3.5e-3 at tol 1e-9. RSS 0.06641 at tol 1e-3. Quiz Q10 (6.25).

## identifiability

- **Corrupted R block repaired** (exercise, C1 solution, "(d) check B3"). Three copies of the `A/k/tb/Jb` lines were there, plus a stray `NA` line and a duplicated `print`. The NA printed as "[1] NA". The block is now the clean 6-line version, matching the Python and Rust twins.
- Fixed in exercise: blocks C2 and C5 call `np.set_printoptions(suppress=True)`. Run in page order, that changed the later C3(a) and C6 output to fixed notation ("[11.89 3.483 1.048 0.]"). The two affected prints now pass `suppress_small=False`, so they print the shown scientific form whatever ran before.
- Fixed in exercise: C2 said `hess_inv` gives SEs "3 times too small". The output shows factors of 2.0, 3.0 and 5.9, so it now says "2 to 6 times too small".
- Figure script: `make-identifiability-figs.R` now uses the pass-2 names: `bateman(t, ka, ke, V)` with a global `D <- 100`, and `tt` for the sample times. The late subset is `tl`/`yl`, and helper arguments are `ts`/`ys`. It runs, and all six SVGs are byte-identical to before.
- Checked, no change: B1–B6 arithmetic (B3 SE 0.195/0.0250, r = 0.546; B5 profile table and Δ = 21.3/0.21/11.4; B6 r = −0.975). Lecture numbers: Wald (0.33, 1.92), profile (0.68, ∞), plateau 1.07, profile vs Wald within 0.005, SE along v3 0.085 (9 %), H vs 2JᵀJ 0.3 %, prediction table (truth 1.86/2.82/1.87/0.57, SEs 0.04 / 0.51, 0.54, 0.05), sloppy 4.5 decades and the 180× ratio, design 0.40 → 0.10, prior SE 0.35. C4 (σ3 0.21, κ 207, 28.6), C6, C7 and quiz Q3/Q6/Q8.
- Not checked: the numbers in the speaker notes of the "Try it" demo slide (SE 0.16 with the 3 h sample moved to 2 h, and so on). They come from the demo's own per-sample noise and cannot be reproduced from the page code.

## constrained-lp

- Code: all R/Python output matches.
- Fixed in exercise: the B8(c) derivation step read "11 RES + FER − (RES + FER + BIO) = 5 BIO − 2g". That equation is wrong, and FER = 2 BIO − 6 RES does not follow from it. The step now reads "(2g + 11 RES + FER) − 2g = 5 BIO − (RES + FER + BIO), i.e. 12 RES + 2 FER = 4 BIO". The result BIO = (2g + 5 RES)/3 was already correct.
- Checked, no change: B1–B7 (B3 p1 = 27/52, misfit 4.923, KKT gradient (−100.9, −51.7, −100.9); B4 corners and prices 2/1, tie at 15; B5 (1.5, 4.5) with prices 3/0/3; B7 10/3 and price 1/3). B8(a)(b)(d). Lecture: half-life 1.45 → 1.33 h (8 %), SSE 0.0023/0.0038, simplex edge rates (both counting conventions), shadow prices 2/1 by re-solving, C(100, 50) ≈ 1e29, FBA 15/5, RES-KO 6.67, oxygen 4 gives growth 10 and acetate 8, slope 5/6. Explore E1–E5 arithmetic. Quiz Q4, Q9 and Q11.

## Renders

All 12 pages render (see below).

## Open

- Nothing open for these units apart from the demo-note numbers above (identifiability) and the literature items in `plan/factcheck.md`.
