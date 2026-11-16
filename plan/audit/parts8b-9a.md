# Pass 3 audit: parts8b-9a (stiffness, phase-plane, gillespie, optimisation-1d)

Method: every R and Python block of lec (R only), exercise and quiz was extracted and run cumulatively per page and
language (`Rscript --vanilla`, `python3 -I`); each shown output line was compared with the real output. All four
figure scripts ran without errors and every referenced figure exists; pass 2 changed no figure data. All by-hand
solutions and quiz answers were recomputed, and the numbers quoted in the text (slides, notes, exercise prose,
the demo answers where they can be computed) were checked against code. No Rust files were touched.

## stiffness

- Code: all R/Python outputs match, except wall-clock times. Those are marked machine-dependent (lec notes, C3 text). In C6 I added "on our machine (yours will differ)" to the "2.5–4 s" claim.
- The lec blocks are not cumulative in page order: "A biology problem first" uses `enzyme` from the next slide, and "Backward Euler on the enzyme" uses `f`/`jac` from the exercise. This is intended and the notes say so. The blocks were run with those definitions supplied.
- Fixed (lec, "Backward Euler" slide): the shown `rbind(...)` output was hand-rounded and not what R prints. The code is now `round(rbind(...), 3)`, so the output shown is the real one (one extra code line).
- Fixed (lec, "Backward Euler on the enzyme"): the same problem with mixed hand rounding. The code is now `round(Y[c(11, 101, 301), ], 4)` and the real output is shown.
- Fixed (lec notes, dimerisation by hand): the Newton errors are 0.38, 0.049, 0.0010, 5e-7 (the notes said 0.00097, 4e-7).
- Fixed (lec notes, "Error against step size"): explicit Euler at h = 1.5 ms gives an error of 4.8e-4, not 4.6e-4. Also, 1.5 ms is not "its largest stable step"; the limit is 2 ms. The sentence is reworded.
- Fixed (lec, "Robertson: steps and time"): "To 10^11 lsoda takes 1595 steps" holds only at `rtol = 1e-8, atol = 1e-12` (default tolerances give 403). The slide now says "(`rtol = 1e-8`)".
- Fixed (exercise B2(b)): "Both are off by similar amounts" was wrong. Backward Euler is off by −0.27 and Euler by +0.62/+0.41. The real numbers are now stated.
- Checked, correct: B1–B7, Q1–Q11, the eigenvalues and time scales (1201, 0.083, ratio 14 400, 2.32 ms, 12 937 steps), 3.37 Newton iterations/step, MM 61 steps and P(10) = 8.24, Robertson max y2 = 3.65e-5, fast eigenvalue −3393, ode45 at atol 1e-10 still 8.89e-6, the BDF2/BE orders, and the damping 3.5e-4 / 2850. The demo answers E1–E5 that can be computed by hand (R factors, limits, step counts ≈ T|λ|/2) are consistent.
- Progress note (deSolve ode45 quirks, machine-specific times): both are stated on the slides/notes; nothing more to do. The stray Chrome profile at `/big/.../chr1` is outside the allowed area; I left it for the user.

## phase-plane

- Code: all R/Python outputs match.
- Fixed (lec notes, bifurcation outlook): "decreasing b from 0.6 to about 0.4 brings the real part below zero" was wrong (Re λ = +0.047 at b = 0.4). It now says "below about 0.35" (Hopf at b = 0.346, B7(c)).
- Checked, correct: J at U on (−1.745, −0.255), τ/Δ, the eigenvalues/eigenvectors (0.934, ±0.357), det −0.4186, saddle (1.094, 1.320) and ±0.1767, Newton biased basins 21.8/57/21.2 % and the ODE's 58/42 % (fig script log), escape times ln(1e4)/0.19, ln(1e7)/0.19 ≈ 85, Sel'kov τ = 0.196, Δ = 0.44, 2π/0.656 = 9.58, period 9.6226, the α-table (u(100) = 1.3717 / 1.1686 against 1.3702 / 1.1051), B1–B8, Q1–Q11, the E4 table (V on / saddle eigenvalues, saddle-node between αv = 2.705 and 2.70, U-on state (2.718, 0.322) at 2.70) and the E5 Re λ / Im λ table.
- Progress note (Sel'kov vs repressilator): reconciled in pass 2. Densest slides: the brief says that no slide overflows. Install: rootSolve is on `intro/install.qmd`.

## gillespie

- Code: all R/Python outputs match (the seeded results reproduce exactly). The EM block shows the NaN warning as a comment; R prints "There were 50 or more warnings", which is acceptable.
- Fixed (lec notes, "Many cells"): "The largest deviation is 1.7 SE (at t = 10)". t = 0.5 deviates as much (1.75 vs 1.70 SE). The notes now say "1.7 SE at t = 0.5 and at t = 10".
- Checked, correct: the by-hand SSA step (0.0602 h, 0.0055 h), the event times 0, 0.133, 0.183, 0.317, the per-event mean 10.5 = (100 + 110)/20, the telegraph g* = 0.1, burst size 22, the Fano factors 16 and 2.67, and the bursty protein Fano 32.47 (solved from the Lyapunov moment equation). Also the tau-leap table (it matches C4's R output), EM NaN 7.6 %, extinction 0.335 ≈ 1/3, the cost table (a0 = 4200, 126 000 steps), B1–B7, Q1–Q11, the E3/E4 exact Fano values and the bursty mean at T = 1 (5.42).
- Progress note (gene parameters vs euler): reconciled in pass 2. Screenshots: the brief says that no slide overflows.

## optimisation-1d

- Code: all R/Python outputs match.
- Fixed (lec notes): log-linear regression gives k = 0.348 only through the origin. The notes now say `lm(log(y) ~ 0 + t)` and give 0.335 with an intercept.
- Fixed (lec notes, 2×): the SSE limit as k → ∞ is 1.41, not 1.47.
- Fixed (lec slide + notes): golden section needs 39 passes (not 38) to get the width below 1e-8, i.e. 41 evaluations, as C1 says.
- Fixed (lec "How fast?" slide + notes): fig05 shows Brent stalling at 4.9e-11, not "near 1e-9 / about 5e-10". The text now says 1e-10 / 5e-11, which agrees with exercise E2 (4.9e-11).
- Fixed (lec "Grid first, then refine"): "121 + 10 evaluations" is now 121 + 9 (counted).
- Fixed (exercise B2(c), B3(b)): the bracket [1.7082, 2.8541] holds after **three** evaluations, not four; the fourth evaluation (2.4164) narrows it to width 0.708.
- Fixed (quiz Q6): `optimize(f, c(0, 10))` does not return 0.3754002 (it returns 0.3753847). The interval is now `c(0, 1)`, which reproduces the stated value.
- Checked, correct: the golden table by hand (0.00216, 0.09150, 0.09417, 0.02204), convexity limit k < 0.80, SSE''(0.7) = 0.25, SSE''(1) = −0.19, the parabolic runs (errors 0.020 … 2e-9; 0.70, −1.51, …, 34), the circadian RSS values (21.46 / 28.54 / 32.92 → 10.36 / 10.31 / 10.25, flat 10.72), the B6 table, the ZTP nll'' = 56.92 and the flat-bottom width 2.0e-8, B1–B6, Q1–Q11.
- No progress-note items for this unit besides the overflow probe (the brief says that no slide overflows).

## Renders

All 12 pages render. Several renders first failed on `.quarto`/site_libs races with other workers; each of these succeeded on a retry. The decks still have 26/26/27/28 slides; none was added. A headless overflow probe after the edits: stiffness (the two edited code slides included) and optimisation-1d have no slide over 720 px. The phase-plane and gillespie decks were changed only in the notes.
