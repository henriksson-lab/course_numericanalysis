# Pass 2 report: parts1-2 (units 1–7)

All title mentions on visible slide/page text now link to `../<slug>/lec.qmd`
(first mention per slide/page and every "Where this leads" bullet); speaker
notes keep plain titles. Part-number-only pointers ("unit 3", "Part 7") on
visible text were replaced by the linked unit title. Exercise R/Python tabs
that call functions defined elsewhere on the page got a one-line
`# uses f() from the Background / the solution of C1` comment (newton already
had these). No code output changed, so no re-runs were needed (except R spot
checks below). Part 1 has no xref files; checks listed here instead.

## intro
- Silent-failure slide: "samples a few points" → "samples 21 points", link to
  *Numerical integration*; notes now state R works up to an upper limit of
  ~4300 and fails from ~4400 (21-point first pass only to ~2300; one bisection
  stretches it), consistent with quadrature/exercise E5 and B5(c). Checked in R
  (2300, 4300 → 0.5; 4400 → 5.7e-5). Exercise C5 solution says the same and
  links quadrature's Explore. (There was no "fine up to ~2300" wording in intro
  before; added the correct statement.)
- Black-box table: Part column now links (3 → lu, 2 → newton (Brent, uniroot
  tol), 9 → gradient-methods, 7 → quadrature, 8 → runge-kutta).
- Instability slide: Part 8 → links to *Initial value problems and Euler's
  method* (which reuses k = 0.8, h < 2.5 h exactly) and *Stiffness*.
- Growth-curve pointers: gradient-methods (not NLS) is the unit that returns
  to this fit and scaling; NLS has multistart. Wording fixed on slides,
  notes and exercise (C4(e), C4(c), Where this leads).
- "Next in Part 1": matches units; added each unit's second biology example
  (dose near the plateau; growth rate from an OD curve). Course-map notes:
  growth curve also returns in convergence-cost and finite-differences.
- Taylor table rows link to newton, error-conditioning, finite-differences,
  euler/runge-kutta, optimisation-1d.

## floating-point
- Links to intro, error-conditioning, convergence-cost, bisection-fixedpoint,
  newton. "variances of large counts" → "variance of genomic positions near
  2.4e8" (error-conditioning's actual CUT&RUN example). Where-this-leads now
  names newton's ZTP/expm1 use.

## error-conditioning
- Links throughout. **Fixed**: "Where this leads" claimed *Newton and secant
  methods* solves a Hill/MM equation for a dose near the plateau as B7's
  ill-conditioned problem; newton doesn't. Now: EC50 by Newton; a flat curve
  makes the step −f/f′ huge (newton's Failure 1), same small slope as B7.

## convergence-cost
- Links. **Fixed**: eigen-algorithms uses 20 000 cells × 2 000 genes, not
  100 000 cells. "Part 2 explains the constant" → newton, |f″/(2f′)| (checked
  = 0.354 for x²−2). Notes: unit slugs replaced by titles.

## bisection-fixedpoint
- Links. **Fixed**: optimisation-1d does not use the harvested-Ricker MSY;
  the bullet now says golden-section would find that maximum and that the unit
  uses its own examples. xref file updated.

## newton
- Links (recap, expm1, Next, Where this leads, quiz Q10). B1(d) notes "4 in
  place of 2". `# Output:` style left as is. xref file updated.

## newton-systems
- **Fixed** recap slide: table now is newton's EC50 run from c0 = 1
  (c_k, f, f′, error 2.26, 0.36, 6.5e-3, 1.6e-6, 9.1e-14 — same as newton's
  "watch the digits" slide; recomputed in R). Slide gained one row; not
  screenshotted — check overflow in audit.
- Structured-sparse has no Jacobians: Where-this-leads reworded. Links
  throughout. xref file updated.

## Needs a change outside this group
- None required. Optional: `quadrature/lec.qmd` notes for slide "Why
  `integrate(dnorm, 0, 20000)` is 0" could mention the 4300/4400 limit (it's
  only in the demo notes and exercise E5); `optimisation-1d` could reuse the
  bisection-fixedpoint MSY (14.26 at N = 46.75) if it wants a back-link.

## Checks
- All 21 pages (lec/exercise/quiz × 7) render (some after retry on
  site_libs races). Decks have 24–27 slide headings (≤ 30). No Rust code
  touched; `cargo build -q --examples` gives no warnings in these files.
