# Cross-references from fem-outlook

- fem-outlook/lec.qmd:24 -> diffusion-explicit, steady-state-2d (Part 10 rectangular grids; About slide)
- fem-outlook/lec.qmd:33 -> diffusion-implicit (Thomas solve, grid refinement as a check; notes)
- fem-outlook/lec.qmd:56 -> diffusion-explicit (mass balance derivation of the diffusion equation; notes)
- fem-outlook/lec.qmd:67 -> splines (piecewise linear interpolation; notes)
- fem-outlook/lec.qmd:81 -> diffusion-explicit (no-flux wall via ghost point; notes)
- fem-outlook/lec.qmd:99 -> finite-differences (second difference (-1, 2, -1)/h^2)
- fem-outlook/lec.qmd:113 -> diffusion-explicit (ghost-point no-flux row; notes)
- fem-outlook/lec.qmd:127 -> reaction-diffusion (length scale sqrt(D/k), synthesis-diffusion-degradation; notes)
- fem-outlook/lec.qmd:174 -> structured-sparse (thomas() function), diffusion-implicit (thomas argument convention; notes)
- fem-outlook/lec.qmd:195 -> diffusion-implicit (grid refinement without an exact solution; notes)
- fem-outlook/lec.qmd:206 -> advection-diffusion (wiggles of central differences at high cell Peclet number; notes)
- fem-outlook/lec.qmd:217 -> steady-state-2d (sparse direct solve), iterative-solvers (3D systems; notes)
- fem-outlook/lec.qmd:229 -> newton-systems (Newton for Navier-Stokes at larger Re; notes)
- fem-outlook/lec.qmd:236 -> advection-diffusion (Poiseuille profile taken as given; T-sensor w = 100 um, mean 1000 um/s)
- fem-outlook/lec.qmd:244 -> advection-diffusion (coupling a computed velocity field into advection-diffusion; notes)
- fem-outlook/lec.qmd:266 -> diffusion-explicit, diffusion-implicit, steady-state-2d, advection-diffusion (Part 10 finite-difference codes)
- fem-outlook/lec.qmd:269 -> nonlinear-least-squares (fitting a model to data, Part 9)
- fem-outlook/lec.qmd:280 -> steady-state-2d (mother-machine exercise: does the model answer the right question; notes)
- fem-outlook/lec.qmd:292 -> iterative-solvers (Next: slide)
- fem-outlook/exercise.qmd:36 -> structured-sparse (Thomas algorithm, Fact 4)
- fem-outlook/exercise.qmd:46 -> structured-sparse (thomas R comment; same in Rust tabs at 144, 366, 529, 692)
- fem-outlook/exercise.qmd:258 -> quadrature (trapezoid rule, Hint 1)
- fem-outlook/exercise.qmd:266 -> diffusion-implicit (thomas argument convention, Hint 3)
- fem-outlook/exercise.qmd:270 -> floating-point (do not test floats with ==, Hint 4)
- fem-outlook/exercise.qmd:790 -> steady-state-2d (Where this leads)
- fem-outlook/exercise.qmd:791 -> iterative-solvers (Where this leads)
- fem-outlook/exercise.qmd:792 -> advection-diffusion (Where this leads: Poiseuille given vs computed)
- fem-outlook/quiz.qmd:111 -> advection-diffusion (Q9: Poiseuille profile used as given)
- fem-outlook/quiz.qmd:118 -> advection-diffusion (Q9 answer)
- fem-outlook/quiz.qmd:130 -> nonlinear-least-squares (Q10: fitting needs many solves, Part 9)

## Running examples defined here

- **Hydrogel slab with uptake**: drug diffusing into a cell-laden hydrogel slab of thickness L = 1000 um on glass, D = 100 um^2/s, uptake k = 0.04 /s, u(0) = 1 (bath), u'(L) = 0; lambda = sqrt(D/k) = 50 um, exact u = cosh((L - x)/lambda)/cosh(L/lambda). 10 linear elements: max nodal error uniform 0.0635, graded x = L (i/N)^2 0.0139. Variant: cell-free top coating of 100 um (k = 0), exact u(100) = 1/(1 + 2 tanh(18)) = 1/3.
- **Pillar chamber**: 600 x 400 um culture chamber with a 3 x 2 array of round pillars, radius 50 um, centres 200 um apart (geometry in make-fem-outlook-figs.R).
