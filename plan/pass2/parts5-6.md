# Pass 2 report: parts5-6 (units 13–17)

All 15 pages (lec/exercise/quiz × 5) render (after retries on glob errors). Deck sizes: power-iteration 25, eigen-algorithms 27, polynomial-interpolation 23, splines 26, fourier-fft 25. `cargo build -q --examples` shows no warnings in these units' files. No Rust code was changed. The only code re-run was `power_rq` on the 3×3 test matrix (R), for the new eigen-algorithms recap slide. All five `plan/xref/<slug>.md` files are tagged.

## power-iteration (unit 13)

Checked: the claims about convergence-cost, bisection-fixedpoint (L/(1−L) false convergence), least-squares-qr, lu (2/3 n³, 2n², never `solve(A) %*% b`), pivoting-conditioning, newton, phase-plane, gillespie and iterative-solvers all hold.

Fixed:
- The Summary "Next" and the exercise "Where this leads" said "PCA of 100 000 cells". eigen-algorithms uses **20 000 cells**, so both now say 20 000, and "huge" became "large". Both are now links.
- Recap-slide notes: the line "high-degree polynomial roots badly conditioned, see How eigen() works" pointed to something that unit does not discuss. The note now says that unit shows what libraries do instead.
- Links: the exercise C5 solution remark and all four "Where this leads" items.
- Solutions: added a one-line pointer saying the R/Python solutions share one session (`L` from the Background, `power_it` from C1, `P` from C4, `C` from C5). Each Rust tab is a complete program.

## eigen-algorithms (unit 14)

Fixed:
- **Recap slide rewritten.** It had been written before power-iteration existed and recapped with its own matrix only. It now recaps power-iteration's actual examples:
  - the fish Leslie matrix (λ1 = 1.177, rate 0.654, 55 steps);
  - the semelparous fish (fails), and the slow-switching cell chain (5406 steps);
  - inverse iteration with a shift, and Rayleigh quotient iteration (cubic for symmetric matrices).

  It keeps the 3×3 `A` as this deck's test matrix, now run with that unit's `power_rq` (51 steps, re-run). The notes say "do not re-teach".
- **RQI is not taught as new.** The "Second trick: shifts" slide and notes now say the Rayleigh shift a_nn *is* the Rayleigh quotient iteration of *Power iteration and friends* (slide and C7 there). Exercise fact 3 says the same, with a link. Slide count is unchanged (27).
- **Wrong back-reference:** "When it stalls" claimed power iteration had failed on a *periodic Markov chain*. In power-iteration the failure is the **semelparous fish**. Fixed in the slide bullet and the notes, and the same wording is used in the Ritz-value notes.
- Links: the About slide, Next, and the four exercise "Where this leads" items.
- Note: there is no URL for the sister course anywhere in the repo, so its mention stays plain text.

## polynomial-interpolation (unit 15)

Checked:
- The OD curve (K = 1.2, N0 = 0.05, r = 0.8/h) and its numbers match convergence-cost and splines: f''(3) = 0.0592, max|f''| = 0.0739 at about 2.27 h, od(2), od(3), od(4).
- The least-squares-qr claims hold: the polynomial fit is a Vandermonde system, kappa grows with the degree, and poly() is a QR of the Vandermonde matrix.
- The fourier-fft claims (trigonometric interpolation, circadian data) and the quadrature claims (Simpson integrates a parabola, Gaussian quadrature) hold.

Fixed:
- Turned first mentions into links on the slides, the exercise page and in "Where this leads": Convergence and cost, Least squares, Splines (x3 in the deck), Fourier, Numerical differentiation, Numerical integration.
- Reworded the "Where this leads" claim that "a finite difference is the derivative of an interpolating polynomial". finite-differences builds its stencils from Taylor series. The new text says that, and adds that each stencil is also the derivative of a local interpolant.
- The "needs ... from C3" code comments (8) now say the functions are defined in the C3 solution on this page. Only comments changed, no outputs.
- Updated plan/xref/polynomial-interpolation.md: all entries tagged.

## splines (unit 16)

Checked:
- The OD parameters match convergence-cost.
- The `thomas(a, b, c, d)` convention matches structured-sparse: all four vectors have length n, and a[1] and c[n] are unused. The convention is now written out on the slide and in Fact 4.
- The convergence-cost V-plot and noise claims, the golden-section claim (optimisation-1d), the ridge claim (svd) and the smoothing claim (fourier-fft) hold.

Fixed:
- Turned first mentions into links on every slide, the exercise page, the quiz and "Where this leads".
- Growth-rate slide: added a linked bullet. *Numerical differentiation* repeats the analysis on its own, denser and less noisy run (every 10 min, sd 0.005), where the smoothing spline has rms 0.0047 and the best wide stencil 0.0053. Those are the numbers that unit actually shows. I shrank the figure from 62 % to 56 % so the slide does not overflow. The notes on slide 2 say the same.
- "Where this leads → Numerical differentiation" said "for functions you can evaluate exactly". It now describes that unit's noisy plate-reader part.
- "Where this leads → Numerical integration" said quadrature integrates a spline; it does not. Reworded: the trapezoid rule is the integral of the piecewise-linear interpolant, Simpson's rule that of local parabolas, and integrating a spline is not treated there.
- Nonlinear least squares and ODE fitting fits a pharmacokinetic (PK) drug model, not the logistic. Reworded the slide notes and "Where this leads". The notes now point to finite-differences' maximum-likelihood fit of the logistic.
- SVD "Where this leads": the svd unit has ridge but never mentions df. Reworded.
- Gene time course: splines uses its own TPM data (eg = 1, 40, 32, ...) at the same times as polynomial-interpolation's C7 fold-change curve g(t), which is a different, smoother curve. Making them agree would mean redoing the figures and the splinelab demo preset. Instead I added one sentence in the recap-slide notes and in exercise C4 saying what differs and why: the overshoot only shows with a steep drop.
- Quiz Q8 now also links Numerical differentiation.
- The "from C2" code comments (R and Python, C3) now point to the C2 solution on the page.
- Updated plan/xref/splines.md: all entries tagged, plus a pass-2 note on the three OD runs.

No code was changed, so nothing needed re-running. Rendered: lec, exercise and quiz for both units (two retries after glob errors, then OK). Decks: polynomial-interpolation 23 slides, splines 26. `cargo build --examples` gives 0 warnings.


## fourier-fft

Checked all 19 xref entries plus the other unit mentions in lec/exercise/quiz. The claims hold:
- the OD curve (K = 1.2, N0 = 0.05, r = 0.8/h) in the fig script and notes matches convergence-cost;
- splines has smoothing splines and a bias–noise slide;
- svd covers ridge and small singular values;
- finite-differences opens Part 7 and shows noise amplification;
- diffusion-explicit has the Fourier-mode solution and the zigzag instability.

Fixes:
- **Links.** The first mention per slide or page is now a link:
  - lec: About slide (Polynomial interpolation, Splines), the Recap slide, the "not periodic" Check (Splines), the image-blur slide (SVD), Next (Numerical differentiation).
  - exercise: all four "Where this leads" items, the C3 solution (least-squares-qr), the C6 solution (Splines).
  - quiz: Q12 answer (SVD).
- **Notes.** The speaker notes now name Convergence and cost (for the OD curve) and Least squares (for the normal equations).
- **"Where this leads" → quadrature.** The quadrature unit shows the trapezoid rule as O(h^2) and never treats the periodic case. I reworded the item so it no longer implies that unit covers it: "the trapezoid rule of that unit is, for smooth periodic functions over a whole period, far more accurate than its usual O(h^2)".
- **Code tabs.** All code tabs were already self-contained. No code changed, so nothing was re-run.
- **Audit notes.** The fourier-fft items (4.8 h, p53, Hes1, the 2-h sampling guideline) are literature fact-checks and are left for pass 3.
- **Checks.** lec, exercise and quiz render. The quiz needed one retry after a glob error. The deck has 25 slides. cargo build gives no warnings in the fourier files.
- The xref file is updated.

Optional outside the group: quadrature could mention the spectral accuracy of the trapezoid rule for periodic integrands, but fourier-fft no longer relies on it.

## Needed outside the group

- **finite-differences/lec.qmd**, slide "The problem: growth rate from a plate reader" (notes): optionally add one sentence saying that *Splines* used a sparser, noisier run of the same curve (15 min, sd 0.01; smoothing-spline rms 0.0067, central differences 0.024), so readers are not surprised by the different numbers. In plan/xref/finite-differences.md the plate-reader data are listed as a "candidate reuse … splines", but splines comes earlier and has its own data; drop that candidate.
- **nonlinear-least-squares**: no change needed. Note for whoever owns it: splines, finite-differences and the splines quiz Q12 describe fitting the logistic as the "model-based alternative", and that unit fits a PK model. It is fine as a general method reference.
- **quadrature** (optional): could mention the spectral accuracy of the trapezoid rule for periodic integrands. fourier-fft no longer relies on it.
- No changes are needed outside the group for power-iteration or eigen-algorithms. Pass 3: the eigen-algorithms recap slide is denser now (two columns), so check it for overflow on screen.
