# Pass 2 report: parts7-8a (units 18–22)

All italic unit titles in `lec.qmd`, `exercise.qmd` and `quiz.qmd` of the five units
are now links (first mention per slide in decks, per page in exercises/quizzes, and
again in "Where this leads"). The xref files are marked `[ok]`/`[fixed: …]`, with a
"Pass 2" section added. Deck sizes have not changed: 26, 27, 26, 26 and 27 slides.
`cargo build --examples` gives no warnings in these files. No Rust code was changed.

## finite-differences
- Richardson form `(4D(h/2) − D(h))/3` is the same as in quadrature's Romberg recap. No change needed.
- [fixed] The lecture's "Smooth first" notes and the exercise's Where this leads said that *Nonlinear least squares and ODE fitting* fits a growth model. It fits a PK model, so the wording now says the method is the same and the example there is a drug-concentration model.
- The recap slide heading now links to *Convergence and cost*. Its content (forward difference, V, √ε step, ratio check) matches that unit.

## quadrature
- [fixed] QUADPACK wording:
  - The lecture's Explore notes said "Click 4400 and 5000: … 2.3e-6". They now say that 4300 works (one bisection), that 4400 and 5000 fail (5.7e-5 and 2.3e-6), and that ~2300 is only the first-pass 5σ rule.
  - The B5(c) solution now adds that the full routine works to ~4300 and fails from 4400 (see E5).
  - Checked in R: 4300 → 0.5, 4400 → 5.69e-5, 5000 → 2.33e-6.
  - No other page says "fine up to 2300".
- [fixed] *Nonlinear least squares and ODE fitting* does not fit chromatography peaks, and it fits its own simulated PK study, not these samples. The peak slide bullet and Where this leads were reworded.
- Checked and correct: the intro C5 example, including the Rust crate missing a peak at 5000; the zero-truncated Poisson MLE 1.2184695 and SE 0.1325 from *Newton and secant methods*; log-sum-exp in *How computers store numbers*.

## monte-carlo
- [fixed] The lecture notes said the OD-curve AUC was "computed with the trapezoid rule in *Numerical integration*". Quadrature used a PK curve, so the wording was changed.
- [fixed] Where this leads said *Identifiability and uncertainty of fitted parameters* teaches the bootstrap. It only mentions it, so the wording was changed.
- Forward reference to *Stochastic simulation* checked and correct: gillespie uses inverse-transform exponential waiting times (`-log(runif)/a0`), a discrete inverse transform via `findInterval(cumsum)` (which matches B2(d)), and Monte Carlo SEs of averages.

## euler
- *Numerical differentiation* (forward difference, rearranged) and *Stiffness and implicit methods* (test equation, Euler bound h|λ| ≤ 2, stability regions, implicit methods) both cover what euler says they do.
- The intro "3. Instability" recap matches, with the same code and output.
- [fixed] The oscillator notes said RK4's larger stability region comes "next unit". runge-kutta does not show stability regions; stiffness does. The notes were reworded.
- [fixed] Where this leads said NLS fits r/K or β/γ. It fits PK rates, so the wording was changed.

## runge-kutta
- SIR parameters agree with euler everywhere: N = 10 000, β 0.3, γ 0.1, (9990, 10, 0), I_max 3007.961.
- [fixed] The Python tabs C2, C3, C5 and C6 are now self-contained: the steppers, `rk_solve` and `adapt_heun` are repeated in each. All four were rerun with `python3 -I` and their outputs are identical to what the page already showed. The R tabs C2, C3 and C6 keep their one-line "uses … from C1/C4" pointer to the same page, which the brief allows.
- [fixed] The recap table uses the error in N(10), while euler's numbers are at t = 5. The notes now say so.
- [fixed] *Equilibria, stability and the phase plane* has no repressilator and no predator–prey model; it uses a toggle switch and Sel'kov. The repressilator notes and Where this leads were reworded.
- lsoda is explained in stiffness (slide "LSODA: switch automatically"). Checked and correct.

## Problems in units outside this group
- **nonlinear-least-squares** (optional): several units (finite-differences, euler) point there for fitting growth or SIR models. NLS uses only the PK model. One sentence there, such as "the same works for a logistic growth curve or SIR", would make those pointers natural.
- **identifiability**: monte-carlo used to point there for the bootstrap. Fixed on our side; no action needed unless identifiability adds a bootstrap.
- **intro** lec notes (slide "Silent failures in two lines") say the rule evaluates "21 points"; that is correct. No change needed. Wording on 2300 vs 4300 does not appear in intro.
