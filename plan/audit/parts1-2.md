# Pass 3 audit: parts1-2 (units 1–7)

Method: every R and Python block from lec (R only), exercise and quiz was
extracted and run cumulatively per page and language (`Rscript --vanilla`,
`python3 -I`, one session per page), and each printed line was compared with
the shown output: inline comments, `# Output:` blocks, and output fences
after the tabset. A reverse check flagged output-like comment lines that no
code produces. All `make-<slug>-figs.R` ran without errors, and every
referenced figure exists. Pass 2 changed no figure data in these units. By-hand
answers, quiz answers and quoted numbers were all recomputed. Decks have
24–27 slides, and no slides were added. No Rust was touched.

Expected differences, not fixed: R warnings, Python RuntimeWarnings and
UserWarnings that the page mentions in a comment ("plus a RuntimeWarning",
"Warning: maxit reached"), and the intended error in the bisection-fixedpoint
quiz. Timings in convergence-cost: the machine ran at load 100–235 during the
audit, so times were 5–150× the shown values. All timing outputs are labelled
as one-machine examples. Text derived from the shown timings (ratios, slopes,
extrapolations) was checked against those shown numbers.

## intro
- Checked: R/Py code (all match), figures, B1–B8, E1–E5, quiz Q1–Q11, slide numbers (doubling times, SSE ratio, log1p bracket, error factors), C1/C4/C5 text.
- Fixed: nothing needed.

## floating-point
- Checked: code (all match), figures, B1–B8, E1–E4, Q1–Q12, slide numbers (470-site underflow with seed, gaps, log-sum-exp, expm1 digits).
- Fixed: exercise C7 text. The error is 11 % at μ = 1e-16 (e^-μ rounds to 1 − 2⁻⁵³) and 100 % only from 1e-17. The old text said 100 % from 1e-16, which contradicted the output.
- Fixed: quiz Q4. The gap ratio is 2¹⁹ ≈ half a million, not "about a million".

## error-conditioning
- Checked: code (all match), figures, B1–B9, E1–E5 formulas, Q1–Q11, slide and notes numbers.
- Fixed: nothing needed.
- Not checked: the demo readouts quoted in E1–E4. They come from the JS demo; only the formulas around them were checked.

## convergence-cost
- Checked: code (all non-timing output matches), figures, B1–B7, E1–E5, Q1–Q12, lecture numbers (h* = 5.4e-8, flop counts, 7.6 h, 17×, memory sizes), C3–C7 text against the shown outputs.
- Fixed: lec "Loops vs vectorised" said "about 17× here" with no remark that the times are machine-dependent (other timing slides have one). Added "(one machine; your times will differ)".
- Fixed: exercise E5 answer listed matrix–matrix and correlation in the wrong order relative to the "×2, 4, 4, 8" and "40 µs, 0.8 s, 11 h, 4.4 h" sentences that follow.

## bisection-fixedpoint
- Checked: code (all match), figures, B1–B8, Q1–Q11, E1–E5 (with the demo's presets), lecture numbers (MSY 14.26, g′ at the roots, Ricker regimes), C1–C7 text.
- Fixed in lec:
  - pole residual "about 10^8" → 10^9 (the code gives −1.07e9);
  - Ricker r = 3 notes: minimum "near 1" → "about 3" (it is 3.2);
  - chaos growth factor "30–100" → "30–120" per 10 years.
- Fixed in exercise:
  - B7: g_b at −1 is not "no convergence". |g′| = 1 is borderline, and this iteration creeps in sublinearly.
  - C6 text: "×1000 per 20 years" → "1000–10 000" (the outputs grow 1000–15 000×).
- Fixed in quiz Q3: x = 1.0000000009 → 0.9999999991. It must lie left of the pole for f ≈ −10⁹, as C1 actually returns.

## newton
- Checked: code (all match), figures, B1–B7, E3–E5, Q1–Q11, lecture tables and failure examples, C1–C6 text.
- Fixed: lec "watch the digits": correct digits "0, 2, 6, 13, 16" → "1, 2, 6, 13, 16" (step 1, 3.6166, already has the leading digit). The matching notes were reworded.
- progress.md note (output convention): newton's trailing `# Output:` blocks were left as they are. They are consistent within the unit and readable, and mass reformatting would only add risk.

## newton-systems
- Checked: code, figures, B1–B8, Q1–Q11, lecture numbers (det J, linear prediction, eigenvalues 0.19 / −2.19, cost 1.5 h per step, basin shares 21/57/21 %).
- Recap slide (the row added in pass 2): recomputed from newton's Background functions from c0 = 1. c_k, f, f′ and errors 2.26, 0.36, 6.5e-3, 1.6e-6, 9.1e-14 are all correct and match newton's "watch the digits" slide. The coordinator checked overflow: none.
- Fixed in lec:
  - "Newton for systems in R" printed with `digits = 16`. In R that gives 2.6180339887498949 0.3819660112501052, not the shown output. Changed to `digits = 15`, which prints exactly the shown line, so "All 16 digits agree" still holds.
  - "Still quadratic" k = 4 error 5.8e-12 → 5.9e-12 (it is 5.864e-12; the notes and exercise already said 5.9e-12).
  - Removed the misleading "(printing switched off)" comment and the "cat line is switched off" note from "Three steady states". The code comment now says "one run per row".
  - Failure 1: (−3.2)^2.5 → (−3.25)^2.5, to match the real step.
- Fixed in exercise C7 (R and Python): the shown table was correct for n = 2, but C6(b) sets the global `n <- 2.5`. Run in page order, C7 gave different roots and maxit counts. Added `n <- 2` / `n = 2` with a comment at the top of C7. After this, the cumulative run matches the shown output exactly. The Rust twin is a separate program and is unaffected.
- progress.md note "Install page: add rootSolve, numDeriv": already present in intro/install.qmd (lines 20–21).

## Render
- All 21 pages rendered (with retries after glob/site_libs races).

## Unresolved
- None. The timings in convergence-cost could not be reproduced under the current machine load, but they are labelled as machine-dependent.
