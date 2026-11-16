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
| 8 | lu | done | running | |
| 9 | pivoting-conditioning | done | done | done |
| 10 | structured-sparse | done | done | done |
| 11 | least-squares-qr | done | done | done |
| 12 | svd | done | done | done |
| 13 | power-iteration | done | done | done |
| 14 | eigen-algorithms | done | done | done |
| 15 | polynomial-interpolation | done | done | done |
| 16 | splines | done | done | done |
| 17 | fourier-fft | done | running | |
| 18 | finite-differences | done | done | done |
| 19 | quadrature | running | | |
| 20 | monte-carlo | done | done | done |
| 21 | euler | done | done | done |
| 22 | runge-kutta | done | done | done |
| 23 | stiffness | running | | |
| 24 | phase-plane | done | done | done |
| 25 | gillespie | running | | |
| 26 | optimisation-1d | done | running | |
| 27 | gradient-methods | running | | |
| 28 | nonlinear-least-squares | running | | |
| 29 | identifiability | running | | |
| 30 | constrained-lp | running | | |
| 31 | diffusion-explicit | running | | |
| 32 | diffusion-implicit | running | | |
| 33 | steady-state-2d | running | | |
| 34 | advection-diffusion | running | | |
| 35 | reaction-diffusion | running | | |
| 36 | fem-outlook | running | (no demo: further reading) | |
| 37 | iterative-solvers | running | (no demo) | |

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
- eigen-algorithms: install irlba, RSpectra. Timings typed into fig script (single run, 1 thread). Matrix 20 000 cells (contentplan says 100 000) — align plan text.
- phase-plane: oscillator = Sel'kov (no repressilator xref existed); runge-kutta may define repressilator — reconcile in pass 2 (fine to have both). Slides not screenshotted; densest separatrix slide, 'Failure: Re λ ≈ 0'. Install: rootSolve.
- runge-kutta: Python C2/C3/C5/C6 tabs import C1/C4 functions from local files (not self-contained) — make self-contained in audit. Repressilator params α=216, α0=0.216, β=5, n=2 unverified vs Elowitz & Leibler. ode_solvers Dopri5 dense-output end-point bug noted.
- fourier-fft: from memory — 4.8-h ultradian example, p53 ~5.5 h, Hes1 ~2 h periods, 'sample every 2 h' guideline. R fft prime-length slowness claim measured here.
