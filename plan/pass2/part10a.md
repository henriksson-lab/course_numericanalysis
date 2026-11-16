# Pass 2 report: part10a (units 31–34)

Group: diffusion-explicit, diffusion-implicit, steady-state-2d, advection-diffusion.

## Grid conventions across the group (summary)

- **diffusion-explicit (31)** defines the reference: nodes x_i = iΔx, i = 0..N, boundary values on the end nodes, µm and s, c scaled 0..1, r = DΔt/Δx², A = tridiag(1, −2, 1), ghost point c_{N+1} = c_{N−1} for a wall. Its lec "About" notes now say that units 32–35 reuse this grid and say where they deviate.
- **diffusion-implicit (32)**: the same grid and the same chamber, now stated. The embryo and PCR-wall models state their differences: N = 200, ghost points at both ends, °C.
- **steady-state-2d (33)**: the same nodes (h = Δx, the same h in x and y), with k = i + (j−1)·nx. It now states two deliberate differences and why:
  - −A, i.e. (D/h²) tridiag(−1, 2, −1), so that the matrix is SPD for Cholesky and CG.
  - Oxygen in µM, because the hypoxia threshold and K_m are absolute.
- **advection-diffusion (34)**: nodes along the channel. Across the channel it uses a cell-centred grid y_j = (j − 1/2)Δy with mirror ghost cells. It now states why: no unknown sits on a wall where u = 0, and the wall is a mirror between cells.
- **Thomas solver**: `thomas(a, b, c, d)` (a = sub-, b = main, c = super-diagonal, a[1] and c[n] unused) is used as stated in structured-sparse, in units 32, 33 and 34.
- **Chamber reuse**: the claims that the chamber is reused were wrong in diffusion-explicit and are now corrected. Only diffusion-implicit reuses the chamber itself. steady-state-2d uses oxygen under a lid, and advection-diffusion uses its own devices.
- **Literature**: tags [G1]–[G4], [T1], [T2], [O1], [O2], [B1], [B2], [H1], [H2], [M1], [M2] all match contentplan.md. Fact-checking against the papers is left for pass 3.
- **Checks**:
  - All 12 pages render (with retries for concurrent glob/site_libs errors).
  - Decks: 25 / 25 / 27–28 / 29 slides, all ≤ 30.
  - `cargo build --examples` gives no warnings in Part 10a files.
  - No code or numbers changed, so nothing was re-run.

## diffusion-explicit (unit 31)

Fixed:
- **FTCS stability link.** The eigenvalue slide said "(Part 5, *Stiffness and implicit methods*)", but stiffness is Part 8. It now reads "eigenvalues from Part 5, the rule from Part 8 ([*Stiffness and implicit methods*])". The notes link Euler's |1 − kh| ≤ 1 ([*Initial value problems and Euler's method*]), the stiffness rule applied eigenvalue by eigenvalue, and `eigen()` ([*How `eigen()` works, and big problems*], Part 5). The About notes link Part 5 via [*Power iteration and friends*].
- **Chamber reuse.** The About notes said "Units 32 to 35 reuse the chamber and the grid". Now they say units 32–35 reuse the grid (and state where they deviate), and only diffusion-implicit reuses the chamber.
- **steady-state-2d description.** The Summary notes said "Steady state in two dimensions takes the chamber to 2D" and the exercise "Where this leads" said "the chamber in 2D". Both were untrue: that unit's example is oxygen under a lid. Both now say it takes the grid to 2D (oxygen under a lid, consumed by cells).
- **Morphogen matrix.** The notes said the method-of-lines matrix is "the same tridiagonal matrix" as structured-sparse's morphogen row. They now say it is the same pattern, which structured-sparse writes with the opposite sign and with degradation k on the diagonal.
- **Citations.** The exercise Background now carries the tags [G1] and [G3] next to the two reviews, as the lecture does.
- **Links.** The first mention on each slide or page, and every "Where this leads" item, are now links in lec, exercise and quiz. Targets: finite-differences, euler, power-iteration, eigen-algorithms, newton-systems, stiffness, fourier-fft, structured-sparse, quadrature, and units 32–35.
- **Verified.** These claims are true and unchanged:
  - finite-differences: the (1, −2, 1)/h² stencil with error h²/12, from the OD growth curve.
  - euler: drug clearance, factor 1 − kh, h = 2 decaying with sign flips, order check on log-log axes, step halving.
  - stiffness: the rule holds for every eigenvalue; the step is set by stability, not accuracy.
  - newton-systems: partial derivatives of F(u, v).
  - fourier-fft: the Nyquist mode.
  - quadrature: the trapezoid rule.
- **Code tabs.** Solution tabs already say "uses … ftcs() and exact() from the Background".
- plan/xref/diffusion-explicit.md is tagged.

## diffusion-implicit (unit 32)

**Fixed**

- Grid made explicit: lec notes of "Recap: FTCS and its limit" and exercise Background "Facts used below" now say it is the node grid of *The diffusion (heat) equation by finite differences* (x_i = iΔx, boundary values on the end nodes, N − 1 interior unknowns, µm and s). Embryo slide notes ("Back to the embryo: a 1D model") and the Background say what differs for the embryo/PCR wall: N = 200 (Δx = 2.5 µm, about 50 points across the transition), no-flux ends via ghost points so the end nodes are unknowns too, temperature in °C. No numerics changed.
- Verified: recap eigenvalue range of A (−π²Δx²/L² … −4, i.e. −0.039 … −40 per s with D/Δx²), FTCS limit r ≤ 1/2, c(100 µm, 60 s) = 0.765, maximum-principle check, compartment reading and ghost points all match diffusion-explicit. `thomas(a, b, c, d)` signature, a[1]/c[m] unused, ≈ 8N operations, diagonal dominance/no pivoting, morphogen row and sqrt(D/k) = 10 cells match structured-sparse. Enzyme + backward Euler with Newton (I − hJ), ode45 riding its stability limit, A-stability, BDF match stiffness. 2D FTCS limit 1/4 is in steady-state-2d notes.
- Richardson extrapolation was attributed to *Numerical integration*; it is introduced in *Numerical differentiation* (quadrature only recaps it). Now links finite-differences, "used again in Numerical integration".
- Exercise "Where this leads": "backward Euler and Crank–Nicolson carry over" to Flow plus diffusion was not true (advection-diffusion uses no Crank–Nicolson); now "backward Euler with one Thomas solve per step carries over, as the step that marches the flow down a channel".
- Links: first mention per slide/notes and all "Where this leads" items now `[*Title*](../<slug>/lec.qmd)` (diffusion-explicit, structured-sparse, stiffness, steady-state-2d, advection-diffusion, reaction-diffusion, quadrature, convergence-cost, finite-differences, fem-outlook, iterative-solvers).
- Audit note (self-contained tabs): every R/Python solution block that calls Background helpers now starts with a comment `# uses thomas, be_step, … from the Background code at the top of this page` (C1–C6).
- Citations [M1], [M2], [H1], [H2] agree with contentplan.md. [H2] is only cited in the exercise Goal (headline-only source; leave to pass 3).
- Deck: 25 slides. plan/xref/diffusion-implicit.md entries tagged.

**Needs a change outside group 31–34:** none found.

## steady-state-2d (unit 33)

Fixed in this unit:
- Grid conventions made explicit. lec About-notes, the "1D first" slide and its notes, and exercise Background fact 2 now each say in one sentence that the grid uses nodes x_i = i h as in diffusion-explicit (h = Δx, the same h in x and y), with unknowns at the interior nodes, a ghost point for a wall, and k = i + (j−1)·nx. They also give the two deliberate differences, with reasons: the matrix is written as minus the second difference, (D/h²) tridiag(−1, 2, −1), so that it is SPD for Cholesky and CG (units 31–32 use tridiag(1, −2, 1)); and oxygen is in µM, not scaled to 0..1, because the threshold and K_m are absolute. The old claim that the conventions were the "same as diffusion-explicit" was only partly true.
- The notes said the 1D matrix was the "same matrix as the morphogen row" in structured-sparse. Corrected to "the same pattern, without the decay term on the diagonal and with fixed values instead of closed ends".
- Audit note on the hypoxia threshold (fig10 slide notes). Removed the claim that the core "changes little" with the threshold. The notes now quote the demo's 5.7 / 10.0 / 16.3 % at 5 / 10 / 20 µM. They add that the density at which a core first appears moves much less, about 1.6–1.7e5 cells/cm² by the constant-uptake estimate (computed by hand, not by running the model). The exercise's E2 solution already said this.
- Exercise "Where this leads" claimed that reaction-diffusion's 2D steps are "a sparse solve like the ones here". They are not, so it now says reaction-diffusion solves the five-point Laplacian with the FFT on a periodic grid.
- Links: the first mention on each slide or page, and every "Where this leads" bullet, now link to the target unit. Code comments stay plain.
- Checked with no change needed: thomas(a, b, c, d) and its a[1], c[n] convention (R and Rust); the structured-sparse 30×30 numbers 4380 / 27 029 / 10 231; CG introduced in structured-sparse; the toggle switch's three roots; lambda = sqrt(D/k) in reaction-diffusion; the 2D FTCS limit 1/4, which is attributed as "derived the same way" and so stays correct; citations [O1], [O2], [B1], [B2], which match contentplan.md ([B1] fact-check left for pass 3).
- Code tabs: R and Python solutions are self-contained or say on the page where a function comes from ("from C2", "from C4"). No code was changed, so no re-runs were needed.
- Checks: lec, exercise and quiz render. The deck has 27 `##` slides plus the title slide. plan/xref/steady-state-2d.md has been updated.

Needs changes outside group 31–34:
- iterative-solvers/lec.qmd, notes on the "Biology: oxygen in a culture chamber" slide (around line 59, "The model is the culture-chamber problem of Steady state in two dimensions"): the problem there is different. It is a 1 mm × 1 mm gel with uniform zero-order uptake (Φ = 12), scaled c with hypoxic below 0.25 (50 µM). steady-state-2d has a 4 × 2 mm chamber under a lid, depth-averaged, at 200 µM, with Michaelis–Menten K_m = 1 µM and hypoxia below 10 µM. Change the sentence to "a variant of the oxygen chamber of Steady state in two dimensions (there 4 × 2 mm under a lid, Michaelis–Menten uptake, hypoxic < 10 µM)". Line 82 ("This is the construction of Steady state in two dimensions") is fine.

## advection-diffusion (unit 34)

Fixed:
- **Grid convention.** One sentence each in the lec notes ("Marching: explicit fails at the walls") and in exercise Fact 6: the cross-channel grid is cell-centred, y_j = (j − 1/2)Δy. This is unlike the node grids of diffusion-explicit and steady-state-2d, so that no unknown sits on a wall where u = 0 (the marching step Δx/u would be infinite there) and the no-flux wall is a mirror halfway between the last cell and its ghost cell (c_0 = c_1). The 1D problems along the channel keep x_i = iΔx. No numerics changed.
- **Mixing length.** I recomputed it with the unit's march()/thomas(), using the criterion "every |c − 1/2| ≤ 0.01" (the demo's 1 % rule):
  - Poiseuille: 3.09 mm at dx = 6.25 µm (the demo's grid), 3.3 mm at dx = 100 µm (the fig10 grid).
  - Plug flow: 4.2 mm, which matches 0.42 ū w²/D.

  The lec notes on the Péclet slide and the H-filter slide used to say "essentially mixed by 4–5 mm, consistent with 4.2 mm". They now give the criterion and flow profile behind each number: plug flow 4.2 mm, Poiseuille about 3.1 mm. The flux fraction is 49.95 % at 5 mm. The biology-slide notes now say "within 1 % by about 3 mm". Exercise Fact 1 now defines "mixed to 1 %" and says the Poiseuille profile mixes earlier (about 0.31 ū w²/D, see E5), so 0.42 ū w²/D is a slightly cautious design estimate. The C4 solution text and quiz Q1 (now marked "plug-flow estimate") were aligned too. E5 already said 3.1 mm.
- **Audit note on E4.** The question now states where C ≤ 2/P comes from: central differences with explicit Euler are stable only if C² ≤ 2r, i.e. u²Δt ≤ 2D, which is C ≤ 2/P because r = C/P. This matches the answer's "negative effective diffusion" explanation.
- **Cross-references.** All verified:
  - diffusion-explicit: FTCS, r ≤ 1/2, the averaging argument, the D table, the checklist, and the zigzag starting at the steep front.
  - diffusion-implicit: backward Euler row −r, 1 + 2r, −r.
  - thomas(a, b, c, d): same argument order in diffusion-implicit and structured-sparse.
  - finite-differences: central difference.
  - optimisation-1d: optimize, minimize_scalar and BrentOpt.
  - identifiability: parameters that appear only as a product.
  - reaction-diffusion: Bicoid and Fgf8.
  - fem-outlook: Poiseuille given here, T-sensor w = 100 µm, 1000 µm/s.
- **Links.** First mentions in the lec, exercise and quiz, and every item in "Where this leads", are now links. The lec notes now also say that thomas comes from structured-sparse. Exercise C3(a) points Python users to solve_banded.
- **Citations.** [G2], [G4], [T1] and [T2] match contentplan.md. Kamholz & Yager 2001 is mentioned in notes under [T1] and hedged; leave it for pass 3.
- **Code tabs and plan files.** The code tabs were already self-contained ("thomas() and march() as in C3" points to the same page). The deck has 29 slides. lec, exercise and quiz render; rendering needed retries because of concurrent glob errors. No code changed, so nothing was re-run. plan/xref/advection-diffusion.md is updated.

Outside the group: none needed. fem-outlook's claims about this unit are correct: Poiseuille taken as given, 1.5 × mean, T-sensor 100 µm / 1000 µm/s, wiggles at high cell Péclet.

## Changes needed outside the group

1. **iterative-solvers/lec.qmd**, notes of "Biology: oxygen in a culture chamber" (about line 59). The notes say "The model is the culture-chamber problem of Steady state in two dimensions", but it is a different problem:
   - iterative-solvers: a 1 × 1 mm gel, zero-order uptake Φ = 12, scaled c, hypoxic below 0.25.
   - steady-state-2d: 4 × 2 mm under a lid, 200 µM, Michaelis–Menten K_m = 1 µM, hypoxic below 10 µM.

   Reword to "a variant of the oxygen chamber of *Steady state in two dimensions*" and name the differences.
2. **reaction-diffusion** (audit note, not checked in depth here): time is in minutes there, against seconds in diffusion-explicit. That unit should state its unit of time and the reason for it in one sentence, as units 32–34 now do for their grids.
