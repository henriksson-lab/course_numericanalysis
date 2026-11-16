# Pass 2 report: part10b-11 (units 35–37)

All cross-references verified against the target units; first mentions per slide/page
turned into links (lec notes included). xref files updated (`plan/xref/<slug>.md`).
No code changed, so no re-runs or figure regeneration were needed. Decks: no slides added.

## reaction-diffusion (35)
- Verified: structured-sparse morphogen row (exchange D = 1/min, k = 0.01/min, source cell 1,
  decay length 10 cells, x[11]/x[1] = 0.368) and `thomas(a, b, c, d)` (same code, same convention).
- Minutes vs seconds: one sentence on slide "Adding a reaction term" (visible), in About notes
  (100 µm²/min = 1.7 µm²/s), and in exercise Background fact 1.
- 2D Turing: notes of "A 2D simulation: IMEX with the FFT" now say why periodic FFT, not the sparse
  Kronecker Laplacian of steady-state-2d (a patch of tissue without edges; the FFT turns the implicit
  solve into a division); exercise fact 7 says the same in a parenthesis.
- "What you can now do" rewritten to match the sidebar: SVD moved to Part 4 group; added Monte Carlo,
  stochastic simulation, constraints, 2D steady states and flow.
- Notes fixes: diffusion-implicit does not use lo/di/up names (reworded to "same thomas call");
  temperature step inside the embryo is ~130 µm wide in diffusion-implicit, not "sharp" (reworded).

## fem-outlook (36)
- Verified: steady-state-2d sparse direct solve and mother-machine exercise (it is exercise C3,
  "Does the model answer the right question?" — now linked to that exercise page); reaction-diffusion
  λ = √(D/k); iterative-solvers for 3D; advection-diffusion T-sensor w = 100 µm, mean 1000 µm/s,
  Poiseuille given, wiggles at cell Péclet > 2; diffusion-explicit mass balance and ghost point;
  diffusion-implicit thomas convention and grid-refinement slide.
- exercise "Where this leads": "the same steady problem" on a rectangle → "the same kind of steady
  problem (oxygen under a lid there)".
- Self-contained tabs: solution tabs C1–C3 (R, Python) now carry a one-line pointer to the Background
  functions (assemble, fem1d, exact, L, D, k, lam).

## iterative-solvers (37)
- Lec notes of "Biology: oxygen in a culture chamber" (former ~l. 59): no longer called the chamber
  of steady-state-2d; now says it is a different problem (1 × 1 mm cell-laden gel, zero-order uptake,
  Φ = 12, scaled matrix) and relates it to steady-state-2d's lid chamber (4 × 2 mm, µM, Michaelis–Menten
  + Newton; same 50 amol/s per cell, same five-point matrix). Exercise Background says the same in one line.
- Kronecker slide notes: "the Kronecker construction of steady-state-2d, on a square and without D/h²".
- Recap of structured-sparse glimpse verified (Jacobi 4159/4158 ≈ 4160, CG 237 ≈ 240, ρ = 0.995,
  [1 2; 2 1] diverges, C8 jacobi()/cg()). Exercise intro now states C8 stopped on an absolute residual
  1e-10 while this unit uses relative 1e-6, and that solutions define their own functions.
- Summary notes / exercise "Where this leads": reaction-diffusion does not produce sparse matrices
  (periodic FFT) — reworded "would, with walls"; mentions steady-state-2d's demo uses SOR.

## Requests for units outside this group
- None required. Optional: steady-state-2d could mention that iterative-solvers uses a different
  (scaled, gel) oxygen problem, if its "Part 11 has the details" note (demo slide) should be precise.
