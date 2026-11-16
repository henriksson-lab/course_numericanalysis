# Cross-references from advection-diffusion

- advection-diffusion/lec.qmd:34 (notes, About this lecture) -> diffusion-explicit (grid, FTCS), diffusion-implicit (backward Euler, Thomas algorithm), optimisation-1d (Brent's method, Part 9), fem-outlook (computing flow fields, Navier-Stokes)
- advection-diffusion/lec.qmd:52 (slide "Recap: diffusion on a grid") -> diffusion-explicit (FTCS update, r <= 1/2, positivity/averaging argument)
- advection-diffusion/lec.qmd:59 (slide "Recap: diffusion on a grid") -> diffusion-implicit (backward Euler row -r, 1+2r, -r; `thomas(a, b, c, d)` with the same argument order)
- advection-diffusion/lec.qmd:70 (notes, Recap) -> diffusion-explicit, diffusion-implicit (positivity argument; Thomas argument order)
- advection-diffusion/lec.qmd:87 (notes, Flow carries the molecules) -> fem-outlook (Stokes/Navier-Stokes for complex geometries, "outlook material")
- advection-diffusion/lec.qmd:106 (notes, Péclet number) -> diffusion-explicit (standard D values table)
- advection-diffusion/lec.qmd:123 (slide "Central differences for dc/dx") -> finite-differences (O(dx^2) central first-derivative difference)
- advection-diffusion/lec.qmd:278 (notes, pulse) -> diffusion-explicit (vectorised update uses old values)
- advection-diffusion/lec.qmd:303 (notes, CFL) -> diffusion-explicit (zigzag appears first where the profile is steep)
- advection-diffusion/lec.qmd:318 (notes, T-sensor 2D) -> steady-state-2d (full elliptic steady problem, sparse 2D solve)
- advection-diffusion/lec.qmd:349 (notes, Marching: explicit fails) -> stiffness (the slowest fluid sets the step: "stiffness again")
- advection-diffusion/lec.qmd:354 (slide "Implicit marching in R") -> diffusion-implicit (`thomas()`)
- advection-diffusion/lec.qmd:375 (notes, Implicit marching) -> diffusion-implicit (r of the implicit diffusion lecture)
- advection-diffusion/lec.qmd:435 (slide "Fitting D") -> optimisation-1d (Brent's method, `optimize`)
- advection-diffusion/lec.qmd:438 (notes, Fitting D) -> identifiability (only D / u-bar identifiable)
- advection-diffusion/lec.qmd:452 (notes, checklist) -> diffusion-explicit (checklist of that unit)
- advection-diffusion/lec.qmd:464 (Summary, "Next") -> reaction-diffusion
- advection-diffusion/lec.qmd:467 (notes, Summary) -> reaction-diffusion (Bicoid, Fgf8), fem-outlook (where flow fields come from)
- advection-diffusion/exercise.qmd:44 (Background) -> diffusion-implicit (Thomas algorithm in the Rust version)
- advection-diffusion/exercise.qmd:281 (Solution B6) -> identifiability (D and u-bar not separately identifiable)
- advection-diffusion/exercise.qmd:693 (C3) -> diffusion-implicit (`thomas()`)
- advection-diffusion/exercise.qmd:1230 (C6) -> optimisation-1d (Brent's method; `optimize` / `minimize_scalar` / `BrentOpt`)
- advection-diffusion/exercise.qmd:1383 (Where this leads) -> reaction-diffusion
- advection-diffusion/exercise.qmd:1384 (Where this leads) -> fem-outlook (computing flow fields, curved walls)
- advection-diffusion/exercise.qmd:1385 (Where this leads) -> nonlinear-least-squares (several parameters), identifiability (B6)
- advection-diffusion/quiz.qmd:82 (Q6 answer) -> steady-state-2d (full steady problem as one sparse linear system)
- advection-diffusion/quiz.qmd:94 (Q7 answer) -> stiffness ("stiffness again")

## Grid conventions (for reconciliation with units 31–33)

Read `diffusion-explicit/lec.qmd` and `diffusion-implicit/make-diffusion-implicit-figs.R` (no lec.qmd or plan/xref
files for units 32–33 existed when this unit was written). Adopted from them: units µm and s; nodes
$x_i = i\Delta x$, $c_i^n$; $r = D\Delta t/\Delta x^2$; `thomas(a, b, c, d)` with a = sub-, b = main, c = super-diagonal,
d = right-hand side (as in the diffusion-implicit figure script); Okabe-Ito colours. New here: the cross-channel
grid for marching is **cell-centred** ($y_j = (j - 1/2)\Delta y$, mirror ghost cells for no flux), so that no
unknown sits on a wall where $u = 0$. If steady-state-2d uses a node-based grid with ghost nodes, say so
when linking.

## Running examples defined here

- **Back-diffusion test problem**: inlet channel $L = 500$ µm, flow $u = 100$ µm/s towards a chamber holding a
  drug ($c(L) = 1$, $c(0) = 0$), $D = 1000$ µm²/s, Pe = 50, layer thickness $D/u = 10$ µm. Exact solution
  $c = \mathrm{expm1}(ux/D)/\mathrm{expm1}(uL/D)$.
- **Stimulus pulse**: 1 s inlet pulse of a protein ($D = 100$ µm²/s) at $u = 200$ µm/s, cells at 2000 µm,
  channel 3000 µm; explicit upwind with $\Delta x = 10$ µm, $C = 0.5$ gives peak 0.6405 (true 0.975).
- **T-sensor / H-filter**: width $w = 100$ µm, $\bar u = 1000$ µm/s, Poiseuille $u = 6\bar u (y/w)(1 - y/w)$,
  sample in $y < 50$ µm; small molecule $D = 1000$, protein $D = 100$, cell $D \approx 0.05$ µm²/s.
- **Christmas tree**: rule $c'_k = (k/n)c_{k-1} + ((n-k)/n)c_k$; 2 or 3 inlets to 8 outlets of 50 µm in a
  400 µm outlet channel, $\bar u = 200$ µm/s, $D = 100$ µm²/s.
- **Synthetic T-sensor profile for fitting**: $D_\text{true} = 400$ µm²/s at $x = 1$ mm, 25 pixels of 4 µm,
  noise sd 0.02; the 25 numbers are printed in exercise C6 and in the figure script.
