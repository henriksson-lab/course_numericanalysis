# Pass 2 report: part9b (units 27–30)

All cross-references in the four units were checked against the target units. First mentions of the titles are now links: per slide in the decks (speaker notes left as plain text), per page and again in "Where this leads" on the exercise and quiz pages. The xref files are updated with `[ok]`, `[fixed: …]` and `[note: …]` tags. Deck sizes: gradient-methods 29, nonlinear-least-squares 26, identifiability 26, constrained-lp 26.

## gradient-methods

- **optimisation-1d does not teach a line search.** It covers f' = 0 and f'' > 0, golden section, Brent, and Newton x − f'/f'' with the check "f must decrease, else shorten the step". It leaves line searches to this unit. Fixes:
  - The slide "What we already have" now recaps the minimum conditions and Newton with its checks.
  - The notes on the About and recap slides, and on the Armijo slide, now say that backtracking turns that check into a rule. They no longer cite a line search from the previous unit.
- **Growth curve vs intro.** Data, start values and results match intro/lec.qmd and exercise C7: NM r = 0.676/min with SSE 2.23, the good fit SSE 0.003, r = 0.524/h. The slide "Back to the growth curve" now also quotes the intro's NM result (r = 0.68, SSE 2.23).
- **Notes, Hessian slide.** The pointer to "the Taylor table" now names intro's "Where Taylor comes back" table.
- **Other targets confirmed:** newton-systems (J convention, `solve(J, -F)`, n+1 evaluations, √ε step, numDeriv, unstable middle state, Broyden), newton (two-site curve, R = 0.5, secant), euler (h < 2/k), pivoting-conditioning, convergence-cost, finite-differences (Richardson).

## nonlinear-least-squares

- All targets were confirmed: least-squares-qr (κ², qr.solve, q_tr_mul), newton-systems (relative FD step), convergence-cost (V-curve), runge-kutta (Dopri5 continuous output), error-conditioning (cancellation), stiffness (lsoda), gradient-methods (LM limits, trust region), constrained-lp (bounds).
- No content changes; links only.

## identifiability

- **Names aligned with nonlinear-least-squares** in lec, exercise and quiz, in R, Python and Rust (9 rust/examples/identifiability_*.rs updated):
  - Model `pk(t, ka, ke, V, D = 100)` → `bateman(t, ka, ke, V)` with global `D <- 100` (Rust: `const D`).
  - Sample times `t` → `tt` (Rust `T` → `TT`).
  - Local helper arguments `tt, yy` → `ts, ys`, so they do not clash with the global names.
  - Re-run: all R and Python blocks, plus the lecture and quiz code, give byte-identical output to before. All 9 Rust examples build without warnings and print exactly the documented output.
  - Added one sentence to the exercise intro naming `bateman`, `tt`, `y`, `D`. The recap slide says the earlier unit's RSS is our SSE.
- **Flip-flop recap** (lecture slide and quiz Q3) now matches nonlinear-least-squares: (1.0151, 0.1906, 20.22) and (0.1906, 1.0151, 3.796), SSE 0.06646744.
- **uniroot.** The notes credited uniroot to bisection-fixedpoint, which only previews it. They now point to newton (where uniroot is introduced), with bisection for the sign-change idea.
- Not changed: make-identifiability-figs.R still uses `pk`/`t`. It already has global `ts`/`ys`, the figures are unchanged, and it is teacher-only.

## constrained-lp

- **Simplex edge rates.** Step 4 of the slide "The simplex idea" now gives both units:
  - carbon edge: −3 per unit x₂ = −1 per unit of proteome freed;
  - back along the proteome edge: −1.5 per unit x₂ = −2 per unit of carbon freed.

  The slide notes, the demo notes and the Explore solution E1 all explain the two ways of counting.
- **"Same null space as in LU" was false.** lu does not treat null spaces. The slide "The stoichiometric matrix" and solution C5 now say that elimination (as in lu) leaves 3 free fluxes.
- **Bounds and standard errors.** The slide "At a bound, the gradient is not zero" claimed that identifiability shows SEs are meaningless at a bound; it does not discuss bounds. The slide now says identifiability's SEs assume an interior minimum.
- **mRNA decay example.** The notes now say that these decay data (A e^{-kt} + c) are not the one-parameter mRNA decay data of optimisation-1d.
- Other targets confirmed: lu (diet 0.55/−0.20/0.65), svd (CD8 −0.09, NNLS/KKT), newton-systems, structured-sparse, gradient-methods (BFGS, L-BFGS-B).

## Problems outside this group

1. **optimisation-1d/lec.qmd, notes of the Summary slide (~line 482):** it says the mRNA decay-rate fit "returns with two parameters (initial amount and rate) in Nonlinear least squares and ODE fitting". It does not: that unit uses the oral PK model. The decay data with a background (three parameters) appear only in constrained-lp, and with different numbers. Suggest: "A three-parameter decay fit with a background comes back in Constraints and linear programming", or drop the sentence.
2. **optimisation-1d/exercise.qmd:297 and :1089** say that line searches come in gradient-methods. That is correct and needs no change; noted for consistency only.
3. **svd/lec.qmd:416** already links constrained-lp for KKT. That is correct.
