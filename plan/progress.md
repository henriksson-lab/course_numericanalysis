# Build progress (teacher-only)

Pipeline per unit: (A) content incl. R/Python/Rust → (B) demo + Explore + slide
→ sidebar entry + glossary merge (coordinator). Then pass 2: cross-references
(plan/xref/*.md → links), pass 3: audit.

| # | Slug | A content | B demo | Sidebar/glossary |
|--|--|--|--|--|
| 1–4 | Part 1 | done | done | done |
| 5 | bisection-fixedpoint | done | running | |
| 6 | newton | running | | |
| 7 | newton-systems | running | | |
| 8 | lu | running | | |
| 9 | pivoting-conditioning | done | running | |
| 10 | structured-sparse | running | | |
| 11 | least-squares-qr | done | running | |
| 12 | svd | running | | |
| 13 | power-iteration | running | | |
| 14 | eigen-algorithms | running | | |
| 15 | polynomial-interpolation | running | | |
| 16–37 | | queued | | |

## Audit notes (pass 3)

- pivoting-conditioning: Cytek "complexity index" = condition number of the reference-spectra matrix — stated from memory, verify or soften.
- least-squares-qr: slides not screenshotted (Chrome hung) — check overflow in audit. exercise.qmd generated from scratch template; xref line numbers may drift.
