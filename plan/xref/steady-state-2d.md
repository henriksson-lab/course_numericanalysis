# Cross-references from steady-state-2d

- steady-state-2d/lec.qmd:34 (notes, About) -> diffusion-explicit (second-difference stencil, grid conventions x_i = i h, units µm/s/µM) [fixed: link; grid sentence now states the two deliberate differences (minus second difference for SPD, µM not scaled)]
- steady-state-2d/lec.qmd:34 (notes, About) -> structured-sparse (tridiagonal and sparse solves, fill-in) [ok, linked]
- steady-state-2d/lec.qmd:34 (notes, About) -> newton-systems (Newton with a Jacobian) [ok, linked]
- steady-state-2d/lec.qmd:70 (Recap) -> diffusion-explicit (second difference, boundary values known) [ok, linked]
- steady-state-2d/lec.qmd:77 (Recap) -> structured-sparse (Thomas O(n), Matrix::solve, fill-in in 2D grid) [ok, linked]
- steady-state-2d/lec.qmd:84 (Recap) -> newton-systems (Newton step, wrong root) [ok, linked]
- steady-state-2d/lec.qmd:93 (notes, Recap) -> structured-sparse (30 x 30 grid: nnz(A) 4380, nnz(L) 27 029 natural, 10 231 reordered); newton-systems (toggle switch, three roots) [ok: 4380 / 27 029 / 10 231 and the 30 x 30 grid match structured-sparse lec; toggle has three steady states]
- steady-state-2d/lec.qmd:107 (notes, steady state = BVP) -> diffusion-explicit (FTCS stability limit; 2D limit 1/4); diffusion-implicit (implicit time steps); Part 8 IVPs (euler) [ok, linked: diffusion-explicit derives only the 1D r <= 1/2; text says 'derived the same way', which is accurate]
- steady-state-2d/lec.qmd:121 (notes, 1D by hand) -> structured-sparse (same matrix as the morphogen row, SPD) [fixed: 'same matrix as the morphogen row' -> same pattern without decay term / closed ends; added why the sign differs from units 31-32; linked]
- steady-state-2d/lec.qmd:144 (notes, 1D in R) -> structured-sparse (thomas function) [ok, linked: thomas(a, b, c, d) from structured-sparse]
- steady-state-2d/lec.qmd:175 (notes, exact solution) -> reaction-diffusion (length scale sqrt(D/k)) [ok, linked: reaction-diffusion uses lambda = sqrt(D/k)]
- steady-state-2d/lec.qmd:331 (notes, sparse direct solve) -> structured-sparse (fill-in, orderings); iterative-solvers (CG, preconditioning, Part 11) [ok, linked]
- steady-state-2d/lec.qmd:376 (notes, Newton slide) -> newton-systems (same method with 3081 unknowns) [ok, linked]
- steady-state-2d/lec.qmd:395 (notes, failure slide) -> newton-systems (several roots, toggle switch) [ok, linked]
- steady-state-2d/lec.qmd:450 (Summary, Next) -> advection-diffusion [ok, linked]
- steady-state-2d/lec.qmd:453 (notes, Summary) -> advection-diffusion (non-symmetric matrix); reaction-diffusion; iterative-solvers (Part 11) [ok, linked iterative-solvers + reaction-diffusion (first in notes)]
- steady-state-2d/exercise.qmd:145 (Background, after tabset) -> structured-sparse (CG introduced); iterative-solvers [ok, linked: CG is introduced in structured-sparse (lec + C8)]
- steady-state-2d/exercise.qmd:231 (Hint 1) -> structured-sparse (Thomas algorithm) [ok: Rust thomas convention a[0], c[n-1] unused matches structured-sparse]
- steady-state-2d/exercise.qmd:302 (Solution B8(c)) -> iterative-solvers [ok, plain (later mention on page)]
- steady-state-2d/exercise.qmd:377 (Rust C1 comment) -> structured-sparse (Thomas algorithm, same signature) [ok: code comment, left plain]
- steady-state-2d/exercise.qmd:1234 (Where this leads) -> advection-diffusion [ok, linked]
- steady-state-2d/exercise.qmd:1235 (Where this leads) -> reaction-diffusion [fixed: reaction-diffusion's 2D uses FFT on a periodic grid, not a sparse solve; wording changed; linked]
- steady-state-2d/exercise.qmd:1236 (Where this leads) -> fem-outlook [ok, linked]
- steady-state-2d/exercise.qmd:1237 (Where this leads) -> iterative-solvers [ok, linked]

## Grid conventions (for reconciliation with diffusion-explicit / diffusion-implicit)

- Same as diffusion-explicit/lec.qmd: x_i = i h, i = 0..N, boundary values known, unknowns at interior points; units µm, s, µM (mM for glucose); no-flux wall by ghost point c_{N+1} = c_{N-1} (closed end is an unknown, last row -2, 2).
- Operator written as MINUS the Laplacian: A = (D/h^2) tridiag(-1, 2, -1) (SPD), equations A c + r(c) = b, boundary values in b with a plus sign.
- 2D: same h in x and y; numbering with x running fastest, k = i + (j - 1) nx (1-based; R column-major, numpy order="F", Rust k = i + nx*j 0-based); A = (D/h^2)(I_ny (x) T_nx + T_ny (x) I_nx).
- diffusion-implicit was not finished when this was written (only figures present); its xref did not exist.
- [pass 2] Checked: diffusion-explicit uses nodes x_i = i dx, tridiag(1, -2, 1), c scaled 0..1; diffusion-implicit the same. The 'same as diffusion-explicit' claim was only partly true; lec About-notes, the 1D slide (+ notes) and exercise Background fact 2 now state the grid in one sentence and the two differences with reasons.
- [pass 2] Audit note 'core size changes little with threshold' softened in lec notes (fig10 slide): area 5.7 / 10.0 / 16.3 % at 5 / 10 / 20 uM; onset density moves only ~1.6-1.7e5 (constant-uptake estimate).

## Running examples defined here

- **Oxygen under a lid** (2D, top view): chamber Lx = 4000 um x Ly = 2000 um, gas-impermeable lid 100 um above a cell monolayer, open along all four edges to air-saturated medium c0 = 200 uM. D(O2) = 2000 um^2/s. Consumption 50 amol/s per cell (assumption) -> R0 = 0.5 uM/s per 1e5 cells/cm^2. Michaelis-Menten Km = 1 uM; hypoxic below 10 uM. Grid h = 50 um: 79 x 39 = 3081 unknowns, nnz(A) = 15169.
  Results: constant uptake min c = 86.15 uM at 1e5 cells/cm^2, -27.71 at 2e5; MM: min 87.19 (1e5), 0.564 (2e5, 10.0 % hypoxic), 0.003 (3e5, 25.3 %); threshold ~1.75e5 cells/cm^2. Plain Newton from c0 at 2e5 converges (12 steps) to a spurious root with min c = -27.61; clipped Newton 8 steps.
  Fill-in (Cholesky, natural vs reordered): h = 50: 240 317 / 43 825; h = 25: 1 984 637 / 257 722; h = 12.5: 16 129 277 / 1 475 954.
- **1D chamber cross-section**: L = 2000 um, N = 4: c = (106.25, 75, 106.25) for R0 = 0.5; (12.5, -50, 12.5) for R0 = 1. Linear uptake k = 0.0025 /s, lambda = 894 um.
- **Mother-machine trench** (1D): L = 25 um, glucose D = 600 um^2/s, Vmax = 0.6 mM/s, Km = 0.005 mM (order-of-magnitude assumptions), rich medium 11.1 mM (0.2 %); closed end by ghost point. Constant-uptake depletion threshold c0 = Vmax L^2/(2D) = 0.3125 mM; L_max at 11.1 mM = 149 um. Closed end (MM, N = 100): 10.79 mM at 11.1 mM.
