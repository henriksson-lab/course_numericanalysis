# Build progress (teacher-only)

Pipeline per unit: (A) content incl. R/Python/Rust → (B) demo + Explore + slide
→ sidebar (python3 plan/sidebar.py; demo labels in plan/demo-labels.tsv) + glossary merge (coordinator). Then pass 2: cross-references
(plan/xref/*.md → links), pass 3: audit.

| # | Slug | A content | B demo | Sidebar/glossary |
|--|--|--|--|--|
| 1–4 | Part 1 | done | done | done |
| 5 | bisection-fixedpoint | done | done | done |
| 6 | newton | done | done | done |
| 7 | newton-systems | done | done | done |
| 8 | lu | done | done | done |
| 9 | pivoting-conditioning | done | done | done |
| 10 | structured-sparse | done | done | done |
| 11 | least-squares-qr | done | done | done |
| 12 | svd | done | done | done |
| 13 | power-iteration | done | done | done |
| 14 | eigen-algorithms | done | done | done |
| 15 | polynomial-interpolation | done | done | done |
| 16 | splines | done | done | done |
| 17 | fourier-fft | done | done | done |
| 18 | finite-differences | done | done | done |
| 19 | quadrature | done | done | done |
| 20 | monte-carlo | done | done | done |
| 21 | euler | done | done | done |
| 22 | runge-kutta | done | done | done |
| 23 | stiffness | done | done | done |
| 24 | phase-plane | done | done | done |
| 25 | gillespie | done | done | done |
| 26 | optimisation-1d | done | done | done |
| 27 | gradient-methods | done | done | done |
| 28 | nonlinear-least-squares | done | done | done |
| 29 | identifiability | done | done | done |
| 30 | constrained-lp | done | done | done |
| 31 | diffusion-explicit | done | done | done |
| 32 | diffusion-implicit | done | done | done |
| 33 | steady-state-2d | done | done | done |
| 34 | advection-diffusion | done | done | done |
| 35 | reaction-diffusion | done | done | done |
| 36 | fem-outlook | done | (no demo: further reading) | done |
| 37 | iterative-solvers | done | (no demo) | done |

## Audit notes (pass 3)

- pivoting-conditioning: Cytek "complexity index" = condition number of the reference-spectra matrix — stated from memory, verify or soften.
- least-squares-qr: slides not screenshotted (Chrome hung) — check overflow in audit. exercise.qmd generated from scratch template; xref line numbers may drift.
- newton: R/Python output shown as trailing `# Output:` block, unlike floating-point's interleaved style — decide one convention in audit. Code on slides renders small; check when recording.
- newton-systems: recap slide assumes *Newton and secant methods* uses √2 (x²−2, x0=1); newton actually uses dose–response + ZTP MLE — fix in pass 2. Slides not screenshotted. Install page: add rootSolve, numDeriv.
- lu: C7 Rust tab omits helper functions (file complete) — only tab not showing whole file. Timings noisy (machine load ~80).
- svd: tool names EPIC, quanTIseq, MuSiC (slide) and CIBERSORT (notes) as constrained-LS deconvolution — from memory, verify. Only 4/27 slides screenshotted. scipy 1.15.2 nnls bug noted in C6 (may go stale). Install: nnls (R).
- least-squares-qr demo: drag interaction untested headless.
- structured-sparse: Thomas signature thomas(a,b,c,d) and morphogen-row example to be reused in Part 10. Rust sprs lacks sparse Cholesky/solve — dense/handwritten substitutes, stated.
- power-iteration: Rayleigh quotient iteration included (slide + C7) — check overlap with eigen-algorithms in pass 2. Cell-state chain 'after Gupta et al. 2011' illustrative. Timing claim 17 s vs 0.3 s single run.
- finite-differences: optim default BFGS gradient = central diff with ndeps=1e-3 — verify. Install: numDeriv.
- splines: tradeSeq / Monocle use spline trends — from memory, verify.
- monte-carlo: from memory — RANDU history, R 3.6 sample() change, ChaCha12 as StdRng, Bellman 'curse of dimensionality', bedtools shuffle / GAT / regioneR. Slides 13–18 figures not seen in browser.
- polynomial-interpolation: figure slides not screenshotted (Chrome hung) — check overflow.
- euler: slides not screenshotted; densest 'Why global error is O(h)', 'Euler by hand'. SIR = N 10000, β .3, γ .1, (9990,10,0), 160 d; runge-kutta told to match.
- eigen-algorithms: install irlba, RSpectra. Timings typed into fig script (single run, 1 thread). Matrix 20 000 cells (contentplan updated).
- phase-plane: oscillator = Sel'kov (no repressilator xref existed); runge-kutta may define repressilator — reconcile in pass 2 (fine to have both). Slides not screenshotted; densest separatrix slide, 'Failure: Re λ ≈ 0'. Install: rootSolve.
- runge-kutta: Python C2/C3/C5/C6 tabs import C1/C4 functions from local files (not self-contained) — make self-contained in audit. Repressilator params α=216, α0=0.216, β=5, n=2 unverified vs Elowitz & Leibler. ode_solvers Dopri5 dense-output end-point bug noted.
- fourier-fft: from memory — 4.8-h ultradian example, p53 ~5.5 h, Hes1 ~2 h periods, 'sample every 2 h' guideline. R fft prime-length slowness claim measured here.
- constrained-lp: from memory — iML1515 2712 reactions / 1877 metabolites, Recon3D ~13 000 reactions. Install: lpSolve, numDeriv.
- diffusion-explicit: slides not visually checked ('FTCS by hand', 'When is the gradient established?', 'No-flux boundary').
- lu demo: E5/swap reveals B5(b)(c) answers; E1/E2 show counts asked in B3(b)/B6(a) — decide in audit. lu/exercise.html needs a re-render (stale vs latest lufactor.js layout — page text unchanged).
- gillespie: gene params differ from euler's (10/min, 0.2/min) — noted in slide notes; reconcile wording in pass 2. Adds SIR early-extinction (2 slides). Only 2 slides screenshotted.
- stiffness: deSolve ode45 quirks on Robertson stated on slides; wall times machine-specific. Stray Chrome profile dir at /big/henriksson/temp/claude-tmp/claude-1001/chr1 (outside scratchpad) — clean up when user present.
- nonlinear-least-squares: install minpack.lm. fig06 zigzag is real (FD step vs ODE tol).
- quadrature: from memory — bioequivalence 80–125 %, regulators accept sampled AUC, WinNonlin/PKNCA rules. Romberg recap assumes Richardson (4D(h/2)−D(h))/3 form — check vs finite-differences in pass 2.
- gradient-methods: slides not screenshotted; densest 'What optim, minimize and argmin do', 'When Newton fails'. From memory: optim ndeps=1e-3 central; scipy forward diff step √ε·max(1,|x|); BFGS gtol 1e-5; McKinnon 1998. rust/examples/gradient_methods_shared/defs.rs included via include! (new pattern).
- identifiability: pk(t, ka, ke, V, D=100) vs NLS unit's bateman() naming — reconcile in pass 2. From memory: Raue 2009, Gutenkunst 2007, DAISY, StructuralIdentifiability.jl. Slides not screenshotted.
- constrained-lp demo: edge rate at optimum −2/−1 per freed variable vs lecture slide '−3 per unit x₂' — same edge, different unit; consider aligning slide wording.
- steady-state-2d: [B1] 'mechanics not nutrients' claim from contentplan, not read — verify. Grid conventions: x_i = i h, µm/s/µM, ghost points, k = i + (j−1)nx — reconcile with diffusion-implicit in pass 2.
- diffusion-implicit: R/Python solution tabs rely on Background-defined helpers (thomas, be_step, cn_step, exact) — decide self-contained convention in audit. h = 10 W/(m²K) from memory. Read [M2], [H1] full text (agent fetched before no-network rule); [M1], [H2] headline only. Embryo 1D model is our own.
- gnlm demo: extra LM stop rule → 13 passes vs lecture slide's 15 (stated on page). Data panel y-axis clips good-start curve at 4 mg/L.
- fem-outlook: notes assume content of steady-state-2d, reaction-diffusion, iterative-solvers, advection-diffusion lecture (T-sensor 100 µm, 1000 µm/s) — verify in pass 2. Hydrogel numbers own round values.
- advection-diffusion: cell-centred grid across channel (vs nodes elsewhere?) — reconcile with steady-state-2d. Device numbers own round values; Kamholz & Yager 2001 cited cautiously. Some code blocks 13–14 lines.
- gradient-methods deck now 29 slides (near ceiling). optrace race/click not tested headless. design demo: typecheck error in demos/src/design/math.ts:109 reported by two workers — fix in audit.
- steady-state-2d lecture notes: 'core size changes little with threshold' — demo shows 5.7/10.0/16.3 % at 5/10/20 µM; soften sentence.
- quadrature: QUADPACK actually works to ~4300 and fails from 4400 (one bisection when error cap hit); '~2300' is the first-pass/5σ rule — check intro + quadrature lecture wording matches. Exercise page horizontal scroll from MathJax.
- advection-diffusion: deck now 29 slides. Mixing length: demo 3.1 mm (Poiseuille, 1 % criterion) vs lecture note '4–5 mm'/4.2 mm plug flow — align wording. E4 uses central-explicit stability C ≤ 2/P (not stated in unit).
- reaction-diffusion: time in minutes vs seconds in diffusion-explicit — reconcile. 2D uses periodic FFT diffusion, not Kronecker Laplacian. From memory: Bicoid λ≈100 µm, D FRAP 0.3 / FCS 7 µm²/s; Fgf8 D≈50, λ≈200 µm; Gregor 2007, Abu-Arish 2010, Berezhkovskii 2010 (notes only). Slides not screenshotted.
- Overflow probe (first run, docs as of pass 1): stiffness slide 9 "Try it: the fast eigenvalue sets the step" overflows (741/720 px). 16 decks gave no result (Chrome timeouts under load) — rerun serially in audit: advection-diffusion, convergence-cost, diffusion-explicit, diffusion-implicit, fem-outlook, floating-point, gradient-methods, identifiability, intro, iterative-solvers, least-squares-qr, lu, newton-systems, optimisation-1d, quadrature, reaction-diffusion, runge-kutta.
- iterative-solvers: built its own oxygen-in-gel chamber (Φ = 12) because steady-state-2d didn't exist yet — align with steady-state-2d in pass 2 or state difference.

## Pass 2 (cross-references) — started

Brief: plan/xref-brief.md. Reports: plan/pass2/<group>.md.

| Group | Units | Status |
|--|--|--|
| parts1-2 | 1–7 | done |
| parts3-4 | 8–12 | done |
| parts5-6 | 13–17 | done |
| parts7-8a | 18–22 | done |
| parts8b-9a | 23–26 | done |
| part9b | 27–30 | done |
| part10a | 31–34 | done |
| part10b-11 | 35–37 | done |
- pass2 parts5-6: check eigen-algorithms recap slide (dense two-column) for overflow. Optional: finite-differences note that splines used sparser/noisier run; quadrature note on trapezoid for periodic functions.
- For part10b-11 worker: iterative-solvers lec notes ~line 59 wrongly call its model the culture-chamber problem of steady-state-2d (it's a 1×1 mm gel, Φ=12) — reword; reaction-diffusion should state minutes vs seconds in one sentence.
- Optional (audit): nonlinear-least-squares could add one sentence that the same methods fit growth/SIR models (several units point there). runge-kutta R tabs C2/C3/C6 still rely on C1/C4 functions (stated) — make self-contained in audit for consistency with Python.
- newton-systems recap slide gained a row — check overflow.

## Pass 3 (audit)

Automated (coordinator):
- Full site render: no warnings; anchored links all resolve.
- Rust: 279/314 tabs match exactly; the rest are machine timings, blank-line formatting of nalgebra matrices, or 17 tabs with interleaved output comments (not auto-comparable). Release vs debug last-digit differences (powi) explained — pages use debug as their headers say.
- Slide overflow: all 37 decks checked, no overflow.
- Literature fact-checks collected in plan/factcheck.md — need web, deferred until course lead is present.

Workers (plan/audit-brief.md, reports in plan/audit/): R/Python re-run, figures, by-hand and quiz arithmetic, text-vs-numbers. Groups as in pass 2.
- audit parts5-6 done. Left: re-run fourier-fft/make-fourier-fft-figs.R on an idle machine (fig04 timings); fourier-fft demo numbers not re-checked.
- audit part10b-11 done. Left: re-time reaction-diffusion and iterative-solvers timings on an idle machine (load ~250 during audit).
- audit parts3-4 done. Left: re-measure lu timing figures (fig02/03) on idle machine; least-squares-qr demo drag untested (needs a browser).
