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
| 12 | svd | done | running | |
| 13 | power-iteration | running | | |
| 14 | eigen-algorithms | running | | |
| 15 | polynomial-interpolation | running | | |
| 16 | splines | running | | |
| 17 | fourier-fft | running | | |
| 18 | finite-differences | running | | |
| 19 | quadrature | running | | |
| 20 | monte-carlo | running | | |
| 21 | euler | running | | |
| 22 | runge-kutta | running | | |
| 23 | stiffness | running | | |
| 24–37 | | queued | | |

## Audit notes (pass 3)

- pivoting-conditioning: Cytek "complexity index" = condition number of the reference-spectra matrix — stated from memory, verify or soften.
- least-squares-qr: slides not screenshotted (Chrome hung) — check overflow in audit. exercise.qmd generated from scratch template; xref line numbers may drift.
- newton: R/Python output shown as trailing `# Output:` block, unlike floating-point's interleaved style — decide one convention in audit. Code on slides renders small; check when recording.
- newton-systems: recap slide assumes *Newton and secant methods* uses √2 (x²−2, x0=1); newton actually uses dose–response + ZTP MLE — fix in pass 2. Slides not screenshotted. Install page: add rootSolve, numDeriv.
- lu: C7 Rust tab omits helper functions (file complete) — only tab not showing whole file. Timings noisy (machine load ~80).
- svd: tool names EPIC, quanTIseq, MuSiC (slide) and CIBERSORT (notes) as constrained-LS deconvolution — from memory, verify. Only 4/27 slides screenshotted. scipy 1.15.2 nnls bug noted in C6 (may go stale). Install: nnls (R).
- least-squares-qr demo: drag interaction untested headless.
- structured-sparse: Thomas signature thomas(a,b,c,d) and morphogen-row example to be reused in Part 10. Rust sprs lacks sparse Cholesky/solve — dense/handwritten substitutes, stated.
