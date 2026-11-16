# Pass 3 audit: part10b-11 (reaction-diffusion, fem-outlook, iterative-solvers)

Method: every R/Python block extracted in page order and run cumulatively per page and
language (`Rscript --vanilla`, `python3 -I`); printed output compared with the shown
comments automatically (whitespace-normalised). Lecture slides that use functions defined
elsewhere (thomas, lcg, ic0, or `assemble` shown before its definition) were run with those
definitions prepended; the splitting fragment (rd lec "Operator splitting") is not runnable
by design. All three figure scripts run without errors; every referenced figure exists.
Machine load was very high (load average ~250) during the audit, so timings were not
reproducible; non-zero counts, iteration counts and all other numbers were.

## reaction-diffusion

- R/Python code: all outputs match except timings (C7: 1.7 s shown, 15-22 s measured under load).
- Fixed: exercise C1 (R and Python) looped `for (N in ...)`, overwriting the Background's
  `N = 50`, so C2/C3 failed when run in page order (subscript out of bounds / broadcast error).
  Loop variable renamed to `n`; outputs unchanged; whole page now runs cumulatively.
- Fixed: C7 text and lecture "2 s in R" now say the timings are machine-dependent.
- Fixed (lec, "When IMEX fails"): "as slow as explicit diffusion on a 30 µm grid" was wrong
  (dx²/(2D) = 2 s at D = 100 µm²/min gives dx = 2.6 µm, a 30 µm grid allows 4.5 min). Now "2.6 µm grid".
- Fixed (quiz Q1): 4.6 λ = 4.605 × 77.46 = 357 µm, not 355.
- By hand B1–B6 and quiz Q2–Q12 re-done: all correct (B2 system and mass balance, B3 G factors,
  B4 commutator, B5 Jacobian/det/d_c = 8.57, B6 q² errors 7.7/2.0/0.5 %).
- Text vs numbers checked: accumulation times 50.2/100.0/150.3, threshold positions
  230.7/138.3/220.2/138.9 (shifts 92 and 10 µm), mass balance 1.32 and 81 %, Strang ratios
  2.7→3.3 (3.7 in fig), fit 200.3 ± 12.2 / 229.9, dispersion 0.272–0.822, q = 0.487, 12.9,
  rate 0.263, spots 47/63/46 and 48/48/47/47, |1 + dt λ| = 1.039: all agree.
- Audit notes (progress.md): minutes vs seconds and periodic-FFT choice already handled in
  pass 2; the rest are literature claims (factcheck.md) or screenshot items, not touched.

## fem-outlook

- R/Python code: lecture and exercise outputs all match (lecture run with `thomas` prepended
  and `assemble` defined before its use on the "Same matrix" slide, which shows output first by design).
- Fixed (lec notes, 2D mesh slide): "as many non-zeros as the node has neighbours (about 7)" →
  "one non-zero per neighbour plus the diagonal (about 7)".
- By hand B1–B3 and quiz Q1–Q10 re-done: all correct (B3 u1 = −0.180, u2 = 0.0611; lumped
  0.0371, 0.0028; √6 λ = 122 µm; Re = 0.1; Poiseuille 1.2 Pa/mm).
- Text vs numbers checked: 72 of 384 grid cells, 2304 triangles, errors 0.0635/0.0139/0.0169
  (ratio 4.6), undershoot 0.072 vs 0.135, N = 4 values −0.17 and lumped 0.037 vs exact 0.0067,
  57.7/50.8 drug amounts, C3 interface errors: all agree.
- Audit note (progress.md): cross-unit assumptions were verified in pass 2; hydrogel numbers
  are own round values. Nothing further.

## iterative-solvers

- R/Python code: exercise outputs all match; lecture outputs match (PCG slide run with the
  exercise's `ic0`: 100/100/40 as in the table; GMRES 99 on the symmetric chamber).
- Fig script (~10 min under load) reproduces every count quoted: Jacobi/GS 705…45 193/22 598,
  CG 24…827, IC 13…275, MG 5, aggregate table (100/203/197/399, 100/203/103/208, 40/73/49/97),
  CG residual 314, GMRES 116 / GMRES(20) 187, centre 0.116, 14 % hypoxic, nnz(L) 27 M / 38 M.
- Fixed (lec notes, "small step" slide): geometric factor ρ/(1−ρ) = 829, not 833.
- Fixed: fill-in table times (11 s, 52 s; measured 7.7 s and 470 s under load) — notes now say
  times are single runs and machine-dependent; quiz Q1 premise and answer marked likewise.
- By hand B1–B5 and quiz Q1–Q12 re-done: all correct.
- Text vs numbers checked: 1910 sweeps/digit, prediction 11 463, κ = 1659 / 26 560, CG bound ~300,
  0.47/h², 6.9 M nnz = 83 MB, Pe h = 0.78, C2 29 / 14.5 sweeps per digit, C3 2.5√κ.
- No audit notes for this unit in progress.md.

## Checks

- No Rust touched. Decks: no slides added.
- All nine pages render (several retries needed for concurrent site_libs/glob errors).
- Deck sizes unchanged: 26 / 19 / 20 slides.
