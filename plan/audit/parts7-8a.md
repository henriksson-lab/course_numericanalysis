# Pass 3 audit: parts7-8a (finite-differences, quadrature, monte-carlo, euler, runge-kutta)

Method: every `r` block (lec, exercise, quiz) and `python` block (exercise, quiz) was
extracted in page order and run cumulatively per page and language (`Rscript --vanilla`,
`python3 -I`). Each top-level expression's printed output was compared with the shown
output (interleaved comments or the trailing `# Output:` block). The only differences
left are harmless: a trailing blank line, a section comment after an output, or an
inline comment on an output line. All `make-*-figs.R` scripts ran without errors, and
every figure referenced in lec/exercise exists. No figure data changed in pass 2 for
these units. All by-hand solutions and quiz answers were recomputed. Numbers quoted in
slide text, notes and solution prose were checked against the code output or recomputed.
No Rust code was changed.

## finite-differences
- R/Python: all lec, exercise and quiz blocks match.
- Figures: script runs; 5/5 figures present.
- By hand B1–B8 (Taylor derivative table, best steps, Richardson, noise sds, forward-mode AD table, gradient cost) and quiz Q1–Q11: all correct.
- Text vs numbers: C1–C6 prose checked (best h, slopes, five-point 5.4e-8, Romberg corner 2.7e-13, scipy 11 evaluations, rms 0.0203/0.0053/0.0041, 85 % bug); lecture's optim slide checked (16 gradient evaluations, 0.17 SE, SEs 0.0035/0.0006).
- [fixed] lec notes said the numDeriv-gradient optim run ends with gradient "1e-4". The real value is 2e-4 to 1e-3, so the notes now say "below 1e-3".
- Audit note on optim: checked empirically by recording the points optim evaluates. Without `gr`, BFGS uses central differences at ±1e-3 (absolute, with parscale 1). The slide and Q9 are correct. numDeriv is on the install page.

## quadrature
- R/Python: all match.
- Figures: script runs; 9/9 present.
- By hand B1–B6 and quiz Q1–Q12: all correct. Checked: log-down areas, λz, 6800/128 intervals, Gauss points, dnorm(43.4) ≈ 10^-409, the 2300 rule, the substitutions, and the 39 %/1.4 % peak height/area figures in Q12.
- Text vs numbers: lec (4 %, −1.2 %, 269 points, 9e-10, Gauss 0.37/8.8/7e-12/0.57, Romberg on √x 1e-3 → 4e-4, 160/962) and C1–C8 prose all match. "Uniform Simpson 513 points for 1e-6" is 1.04e-6, which is correct to the one digit shown. Left as is (the figure annotation says the same).
- [fixed] C4(b) solution: the chromatogram value at the nearest Simpson point (40 s = 13σ from the peak) is 2.5e-36, not 1e-35. It now says 3×10^-36, which is consistent with the printed area 1.99e-34.
- [fixed] Audit note "exercise page horizontal scroll from MathJax": the Background display equation also held three AUC values on one line. It now holds only C(t); the AUC values follow in a text line (values re-checked with integrate).
- Literature items (bioequivalence, WinNonlin/PKNCA) left for factcheck.md. The QUADPACK 4300/4400 wording was already fixed in pass 2 and is consistent.

## monte-carlo
- R/Python: all match. The lec RANDU snippet uses `randu()` from exercise C1, as the slide says, and has no shown output. With C1's randu the "only −5..9" claim holds.
- Figures: script runs; 8/8 present.
- By hand B1–B7 and quiz Q1–Q11: all correct, including the seed-42 R and numpy numbers in Q1.
- Text vs numbers: lec (0.0131 vs 0.0122, 3 h / 3 million years, 680 000, seeds 2–6 estimates, SE range 0.011–0.029 for seeds 1–5, ESS 3306/65 997, 1.2 % max weight, 1.364^d, 78 %/0.5 %/2 %, 0.75 %) and C1–C8 prose checked against the R, Python and Rust outputs.
- Literature items (RANDU history, R 3.6 sample(), ChaCha12, Bellman, bedtools/GAT/regioneR) left for factcheck.md.

## euler
- [fixed] lec slide "Euler is unstable for large h (logistic)": the shown output was `1.432 0.770 1.432 0.767`, but R prints `1.431 …`, so it now shows 1.431. Also fixed the spacing of the `log10(h)` output on "Measure the order".
- Python: all match.
- Figures: script runs; 8/8 present.
- By hand B1–B8 and quiz Q1–Q12: all correct. Recomputed: the B3 steps, the B4 factors, the B5 powers, the B7 SIR step, the B8 amplitude identity, the Q9 2-cycle 0.768/1.432, and the Q10 peak +17.5 % ("18 %").
- Text vs numbers: lec (24 %, local-error table, 0.17 h, 8300/830 000 steps, SIR 3127/3020, oscillator 3.4 %/2.2/2850/8.6 %) and C1–C5 prose and E1–E5 answers checked.
- SIR parameters agree with runge-kutta.

## runge-kutta
- [fixed] The R tabs C2, C3 and C6 are now self-contained. C2 and C3 repeat the steppers and `rk_solve` from C1 (C2 also defines `methods`), and C6 repeats `adapt_heun` from C4, with the same "repeated so this block runs on its own" comment the Python tabs use. Every R and Python tab of the exercise was rerun on its own in a fresh process: all outputs are identical to the page.
- R/Python: all lec, exercise and quiz blocks match, both cumulatively and standalone.
- Figures: script runs; 10/10 present.
- By hand B1–B7 and quiz Q1–Q12: all correct. Recomputed: the RK4 table, the 31× ratio, the B3 orders and costs (14 000/430), the B4 controller (0.159, err 0.81), the B6 V values, the B7 interpolated peak 38.88, and lynx above 20 for 0.80 of 5.48 years.
- Text vs numbers: the cost figures 61 000/1 400/140 were re-found (60 992/1 426/140). Also checked: the repressilator range 1.25–59 and period 8.69; RK4 h = 0.006 → 2.9e-4 at 33 000 evaluations; ode45 160/962 and 67 steps; the fig08 atol floors (3.3–5 % and 0.19 %); and the solve_ivp atol = 10 claim (2.8–4.7 %).
- [fixed] C3 solution prose: "tolerance 1e-3 drifts by 0.02–0.04". The real values are 0.0167 (R), 0.0165 (Python) and 0.0387 (Rust), so it now says 0.017–0.04.
- Repressilator parameters and the Dopri5 note are literature/crate items, left for factcheck.md.

## Renders
All 15 pages (lec, exercise, quiz × 5 units) were rendered with `quarto render`. Deck sizes are unchanged (no slides added).

## Not resolved / for the coordinator
- Explore-answer numbers that come only from the JS demos (noisyderiv, eulerstep, rk demo, mclab, quadrature demo) were spot-checked where a closed form exists (E3 noise formula, eulerstep E3/E5). They were not re-run in a browser.
- The remaining literature items for these units are listed in factcheck.md.
