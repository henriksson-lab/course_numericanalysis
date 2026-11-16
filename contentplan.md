# Content plan — Numerical analysis for biologists

A working draft of the course content. This is a planning document, not a
contract: units shift, merge and split as the
[method ↔ biology pairings](designprinciples.md) become clearer. See
[audience.md](audience.md) for who we are writing for.

Open questions are marked **[?]**.

## Scope: what the standard textbooks cover

The topic list was checked against the tables of contents of five common
textbooks (from their published tables of contents, not a page-by-page read):

| | Topic | Burden & Faires *Numerical Analysis* | Sauer *Numerical Analysis* | Heath *Scientific Computing* | Driscoll & Braun *Fundamentals of Numerical Computation* | Trefethen & Bau *Numerical Linear Algebra* |
|--|-------|:--:|:--:|:--:|:--:|:--:|
| A | Floating point, error, conditioning, stability | 1 | 0 | 1 | 1 | III |
| B | Nonlinear equations (1 variable) | 2 | 1 | 5 | 4 | — |
| C | Linear systems: elimination, LU, pivoting | 6 | 2 | 2 | 2 | IV |
| D | Iterative linear solvers (Jacobi, CG, GMRES) | 7 | 2, 4 | 11 | 8 | VI |
| E | Least squares, QR, SVD | 8 | 4, 12 | 3 | 3, 7 | I–II |
| F | Eigenvalues: power iteration, QR algorithm | 9 | 12 | 4 | 7, 8 | V |
| G | Nonlinear systems (Newton in $n$ dimensions) | 10 | 2 | 5 | 4 | — |
| H | Interpolation, splines | 3 | 3 | 7 | 5, 9 | — |
| I | Approximation, orthogonal polynomials, Chebyshev | 8 | 3 | 7 | 9 | — |
| J | Fourier, FFT | 8 | 10, 11 | 12 | 9 | — |
| K | Numerical differentiation | 4 | 5 | 8 | 5 | — |
| L | Quadrature | 4 | 5 | 8 | 5 | — |
| M | ODE initial value problems, stiffness | 5 | 6 | 9 | 6 | — |
| N | Boundary value problems | 11 | 7 | 10 | 10 | — |
| O | PDEs (finite differences, FEM) | 12 | 8 | 11 | 11–13 | — |
| P | Optimisation | 10 (partly) | 13 | 6 | 4 (nonlinear LS only) | — |
| Q | Random numbers, Monte Carlo, SDEs | — | 9 | 13 | — | — |

**What we keep, and how heavily**, given the audience (biologists, one
calculus + one linear algebra course, applied focus):

- **Core:** A, B, C, E, F, H, K, L, M, P — every biologist using a numerical
  library meets these.
- **Light (one unit or folded in):** D (only the idea: why iterative solvers
  exist for big sparse problems), G (Newton for systems, as the bridge to
  optimisation), I (Chebyshev only as the cure for Runge), J (one unit).
- **Promoted beyond a standard NA course** because biology needs them:
  Q (Monte Carlo, Gillespie stochastic simulation), parameter fitting of ODE
  models, identifiability, log-space arithmetic for likelihoods.
- **Core, finite differences only:** N and O, restricted to diffusion, heat
  conduction and advection–diffusion on rectangular grids (one part, five
  units), set in **microfluidics**. Boundary value problems enter as the
  steady state of a diffusion problem rather than as a separate topic.
- **Further reading only:** finite elements and complex geometries, flow
  (Navier–Stokes) solvers, iterative solvers in depth, convergence proofs in
  general.

## Delivery

- Self-paced, no time limit. Units in sidebar order; no calendar.
- One unit = one folder = `lec.qmd` (≤ 30 slides) + `quiz.qmd` + `exercise.qmd`,
  usually with one interactive demo.
- **R is the main language.** Slides are R only. Exercise and quiz pages
  show every code example, hint and solution in R, Python and Rust (tabs).
- **[?]** Examination / certificate: is there an exam (written, by hand) or
  hand-in projects? This decides how much by-hand drill goes in the exercises.

## Units

Biology pairing in *italics*. Slugs are the planned folder names.

### Part 1 — Numbers on a computer

1. **`intro` — Why numerical analysis?** What a numerical method is; three
   ways a computation goes wrong (rounding, truncation, instability); error vs
   cost. Recap of Taylor polynomials as the workhorse of the course.
   *A growth-curve fit that returns a confident wrong answer.*
2. **`floating-point` — How computers store numbers.** Binary, IEEE double,
   machine epsilon, rounding, overflow/underflow, `NaN`/`Inf`; why
   `0.1 + 0.2 != 0.3`; comparing floats with a tolerance.
   *Multiplying 10 000 probabilities underflows to 0: log-likelihoods,
   log-sum-exp (as in HMMs and phylogenetics), `log1p`/`expm1` for small rates.*
3. **`error-conditioning` — Error, conditioning and stability.** Absolute vs
   relative error, significant digits, catastrophic cancellation, condition
   number of a problem vs stability of an algorithm.
   *The one-pass variance formula on large counts loses all digits; computing
   fold changes of nearly equal numbers.*
4. **`convergence-cost` — Convergence and cost.** Truncation vs rounding error,
   order of accuracy, reading orders off log–log plots, big-O operation counts,
   timing code.
   *Why a method fine for 100 genes is useless for 20 000 genes × 100 000 cells.*

### Part 2 — Solving nonlinear equations

5. **`bisection-fixedpoint` — Bisection and fixed-point iteration.** Bracketing,
   guaranteed convergence, stopping criteria; fixed-point iteration, cobweb
   plots, the derivative condition $|g'(x^*)| < 1$, linear convergence.
   *Steady states of a discrete-time population model (Beverton–Holt,
   Ricker): when iteration converges, oscillates, or goes chaotic.*
6. **`newton` — Newton and secant methods.** Derivation from Taylor, quadratic
   convergence, secant method, failure modes (bad start, flat derivative,
   cycles), safeguarded methods (Brent: `uniroot`, `brentq`).
   *Solving a Hill / Michaelis–Menten equation for the dose giving 50 %
   response; a maximum-likelihood estimate as the root of the score.*
7. **`newton-systems` — Newton for systems of equations.** Partial derivatives,
   the Jacobian (introduced here), one linear solve per step, finite-difference
   Jacobians.
   *Steady states of a two-gene toggle switch; multiple solutions =
   bistability.*

### Part 3 — Linear systems

8. **`lu` — Gaussian elimination as LU.** Elimination recap, LU factorisation,
   forward/back substitution, factor once and solve many times, $O(n^3)$ cost,
   why never to compute `solve(A) %*% b`.
   *Stoichiometric balance / tracer mixing with many right-hand sides.*
9. **`pivoting-conditioning` — Pivoting, norms and conditioning.** Why
   elimination without pivoting fails, partial pivoting, vector and matrix
   norms, condition number, residual vs error, rule of thumb "lose
   $\log_{10}\kappa$ digits".
   *Spectral unmixing (flow cytometry / fluorescence) with overlapping
   spectra: an ill-conditioned matrix amplifies noise.*
10. **`structured-sparse` — Structure: symmetric, banded, sparse.** Cholesky for
    SPD matrices (covariance), tridiagonal systems in $O(n)$, sparse storage
    (`Matrix`, `scipy.sparse`), fill-in; a glimpse of iterative solvers
    (Jacobi, conjugate gradient) for huge systems.
    *Cell × gene count matrices are > 90 % zeros; kinship / covariance
    matrices are SPD.*

### Part 4 — Least squares

11. **`least-squares-qr` — Least squares: normal equations vs QR.**
    Overdetermined systems, the normal equations, why they square the condition
    number, Householder/Gram–Schmidt QR, what `lm` / `numpy.linalg.lstsq`
    actually do.
    *Calibration (standard) curves for qPCR and protein assays; polynomial
    fits that go badly conditioned.*
12. **`svd` — The SVD, rank deficiency and regularisation.** SVD as rotate–
    stretch–rotate, numerical rank, pseudoinverse, truncated SVD, Tikhonov
    (ridge) regularisation; non-negative least squares.
    *Deconvolving bulk RNA-seq into cell-type proportions; collinear
    predictors.*

### Part 5 — Eigenvalue problems

13. **`power-iteration` — Power iteration and friends.** Power method,
    convergence rate $|\lambda_2/\lambda_1|$, Rayleigh quotient, inverse and
    shifted iteration.
    *Leslie matrix: long-term growth rate and stable age distribution;
    stationary distribution of a Markov chain.*
14. **`eigen-algorithms` — How `eigen()` works, and big problems.** The QR
    algorithm in outline, symmetric vs nonsymmetric, sensitivity of
    eigenvalues; Krylov/Lanczos methods for a few eigenpairs of a large sparse
    matrix (`irlba`, `scipy.sparse.linalg.eigsh`).
    *PCA of a 20 000-cell single-cell dataset without forming the covariance
    matrix.*

### Part 6 — Interpolation and approximation

15. **`polynomial-interpolation` — Polynomial interpolation.** Lagrange and
    Newton forms, interpolation error, Runge phenomenon, Chebyshev nodes;
    interpolation vs fitting.
    *Reading off values between sparse time points; why a high-degree
    polynomial through a growth curve is a bad idea.*
16. **`splines` — Splines.** Piecewise linear, cubic splines, boundary
    conditions, smoothing splines and the bias–noise trade-off.
    *Smoothing a gene's expression along a time course or pseudotime; growth
    rate as the derivative of a smoothed OD curve.*
17. **`fourier-fft` — Fourier series and the FFT.** Trigonometric
    interpolation, DFT, $O(n\log n)$ FFT, sampling, Nyquist, aliasing,
    filtering, convolution theorem.
    *Detecting circadian rhythms in time-series expression; smoothing a noisy
    trace; image blur as convolution.*

### Part 7 — Numerical differentiation and integration

18. **`finite-differences` — Numerical differentiation.** Forward, backward and
    central differences from Taylor, truncation vs rounding trade-off and the
    optimal step, Richardson extrapolation; automatic differentiation in one
    slide (what `torch`/`jax`/Stan do instead).
    *Instantaneous growth rate from a measured curve; gradients for fitting.*
19. **`quadrature` — Numerical integration.** Riemann, trapezoid, Simpson,
    composite rules, error order, Romberg, adaptive quadrature, Gaussian
    quadrature, what `integrate`/`quad` do, infinite intervals and
    singularities.
    *Area under the curve in pharmacokinetics (AUC) and chromatography peaks;
    normalising a probability density.*
20. **`monte-carlo` — Random numbers and Monte Carlo.** Pseudo-random
    generators and seeds, inverse-transform sampling, Monte Carlo integration,
    $1/\sqrt{N}$ error, curse of dimensionality, variance reduction (importance
    sampling).
    *Probability that a randomly placed set of binding sites overlaps a gene;
    a marginal likelihood in many dimensions.*

### Part 8 — Ordinary differential equations

21. **`euler` — Initial value problems and Euler's method.** ODEs as models,
    direction fields, Euler, local vs global error, order 1, instability for
    large steps; systems of ODEs; higher-order ODEs as first-order systems.
    *Exponential and logistic growth; SIR epidemic.*
22. **`runge-kutta` — Runge–Kutta and adaptive step size.** Midpoint/Heun,
    classical RK4, embedded pairs and error control (RK45 / Dormand–Prince),
    tolerances `rtol`/`atol`, dense output, events; `deSolve::ode` and
    `scipy.integrate.solve_ivp`.
    *Lotka–Volterra predator–prey; the repressilator.*
23. **`stiffness` — Stiffness and implicit methods.** Stability region, the
    test equation $y' = \lambda y$, backward Euler, BDF, why implicit methods
    need Newton (ties back to unit 7), LSODA auto-switching.
    *Enzyme kinetics with fast binding and slow catalysis; Robertson's
    chemical kinetics problem.*
24. **`phase-plane` — Equilibria, stability and the phase plane.** Nullclines,
    fixed points, linearisation via the Jacobian, eigenvalues decide stability
    (ties back to Part 5), limit cycles; numerical continuation in one slide.
    *Bistable genetic switch, oscillating cell-cycle / circadian models.*
25. **`gillespie` — Stochastic simulation.** When copy numbers are small,
    chemical master equation in words, Gillespie's SSA, tau-leaping, comparing
    with the ODE mean; Euler–Maruyama for SDEs in brief.
    *Transcriptional bursting and cell-to-cell variability in gene expression.*

### Part 9 — Optimisation and model fitting

26. **`optimisation-1d` — Optimisation in one dimension.** Minimum vs root,
    golden-section search, Brent, Newton for minimisation, convexity, local vs
    global minima.
    *Maximum-likelihood estimate of a single rate parameter.*
27. **`gradient-methods` — Gradient descent, Newton and quasi-Newton.**
    Gradient and Hessian (introduced here), steepest descent and line search,
    Newton, BFGS / L-BFGS, derivative-free Nelder–Mead; scaling and
    starting values; what `optim` / `scipy.optimize.minimize` do. Stochastic
    gradient descent in one slide (as used to train neural networks).
    *Fitting a dose–response (4-parameter logistic) curve.*
28. **`nonlinear-least-squares` — Nonlinear least squares and ODE fitting.**
    Gauss–Newton, Levenberg–Marquardt (`nls`, `minpack.lm`,
    `scipy.optimize.least_squares`), fitting ODE parameters to time-course data,
    log-parameterisation for positive rates, multistart.
    *Fitting a pharmacokinetic or viral-dynamics ODE model to measurements.*
29. **`identifiability` — Identifiability and uncertainty of fitted
    parameters.** Structural vs practical identifiability, the Hessian /
    Fisher information and standard errors, profile likelihood, sloppy models.
    *Two parameters that only appear as a product; why a model fits the data
    but its parameters cannot be trusted.*
30. **`constrained-lp` — Constraints and linear programming.** Bounds,
    Lagrange multipliers in one picture, linear programming and the simplex
    idea, LP solvers.
    *Flux balance analysis of a metabolic network.*

### Part 10 — Diffusion and heat: finite-difference PDEs

The running example is a **microfluidic chip**: channels and chambers a few
hundred µm wide, where flow is laminar, mixing happens by diffusion, and
temperature and concentration fields decide what the cells see. Every unit
reuses earlier parts: Taylor (Part 1) for the difference stencils, sparse and
tridiagonal solves (Part 3) for implicit steps, stiffness (Part 8) for the
time step limit, eigenvalues (Part 5) for the stability analysis.

**No in-house data.** Each unit is anchored in a published device (sources
below). We rebuild a simplified version: geometry and parameters taken from
the paper or from standard values (diffusion coefficient of a small molecule
in water ≈ 10⁻⁹ m²/s, of a protein ≈ 10⁻¹⁰–10⁻¹¹ m²/s, of O₂ ≈ 2·10⁻⁹ m²/s),
and any "measurements" are **synthetic**, generated from the model plus noise.
Say so on the slide. Check every number against the paper before use (the
sources were found by a literature search that only read abstracts).

31. **`diffusion-explicit` — The diffusion (heat) equation by finite
    differences.** Fick's law and the diffusion equation (the same equation as
    heat conduction), discretising space and time on a grid, the central
    second-difference stencil, the explicit FTCS scheme, the stability limit
    $D\,\Delta t/\Delta x^2 \le 1/2$ and what violating it looks like, boundary
    conditions (fixed concentration / temperature vs no-flux), the
    diffusion time scale $L^2/D$.
    *A flow-free chemotaxis chamber [G1, G3]: a channel between a source and
    a sink reservoir. How long until the linear gradient is established, and
    why does doubling the channel length quadruple the wait? Exercise: a
    small molecule vs a protein chemokine.*
32. **`diffusion-implicit` — Implicit schemes and the method of lines.**
    Backward Euler and Crank–Nicolson, one tridiagonal solve per step in
    $O(n)$, accuracy in time vs space, method of lines (the grid turns the PDE
    into a large stiff ODE system, solved with Part 8's tools); checking a
    solver against an exact solution and with grid refinement.
    *Heat: a Drosophila embryo held between two laminar streams at different
    temperatures, so its two halves develop at different speeds [M1, M2]. How
    sharp is the temperature step inside the embryo? Exercise: heat
    conduction through the wall of a PCR chip, and why chip thickness and
    material matter [H1, H2].*
33. **`steady-state-2d` — Steady state in two dimensions.** Steady diffusion
    as a boundary value problem (Poisson / Laplace equation), the 2D five-point
    Laplacian, assembling it as a sparse matrix (Kronecker products), sparse
    direct solve, sources and sinks; nonlinear (Michaelis–Menten) uptake
    solved with Newton (ties back to Part 2).
    *Oxygen in a culture chamber where the cells' own respiration creates a
    hypoxic gradient [O1, O2]: where do the cells go hypoxic, and how does it
    change with cell density? Exercise (1D): nutrient in the dead-end
    trenches of a "mother machine" for bacteria [B1, B2]. Do cells at the
    closed end starve? (The paper finds mechanics, not nutrients, limits
    growth: a good discussion of whether the model answers the right
    question.)*
34. **`advection-diffusion` — Flow plus diffusion.** Advection–diffusion
    equation, the Péclet number, upwind vs central differences, numerical
    diffusion and wiggles, the CFL condition, a fixed (Poiseuille) flow
    profile given rather than solved for.
    *The T-sensor / H-filter [T1, T2]: two streams side by side; small
    molecules cross, cells and proteins do not. Then the "Christmas-tree"
    gradient generator [G2] used to expose neutrophils to IL-8 gradients
    [G4]. Exercise (ties to Part 9): estimate a diffusion coefficient by
    fitting the model to a synthetic downstream profile, as the original
    T-sensor papers did [T1].*
35. **`reaction-diffusion` — Reaction–diffusion and pattern formation.**
    Adding reaction terms, the synthesis–diffusion–degradation model and its
    length scale $\lambda = \sqrt{D/k}$, operator splitting / IMEX (diffusion
    implicit, reaction explicit), 2D simulations, sensitivity to grid and time
    step; Turing patterns as an outlook.
    *The Bicoid morphogen gradient in the Drosophila embryo, and the
    microfluidic temperature-step experiment that perturbs it [M1]; embryo
    trap arrays for imaging the Dorsal gradient in many embryos [M3];
    Fgf8 in zebrafish, shaped by diffusion and receptor-mediated
    degradation [M4].*

#### Literature sources for Part 10

Gradients and chemotaxis
- [G1] Kim, Kim & Jeon (2010). Microfluidic gradient generation for 2D and 3D cell environments (critical review). *Integrative Biology*. <https://pubs.rsc.org/en/Content/ArticleLanding/2010/IB/C0IB00055H>
- [G2] Jeon, Dertinger, Chiu, Choi, Stroock & Whitesides (2000). Generation of solution and surface gradients using microfluidic systems. *Langmuir* 16:8311–8316.
- [G3] Review of microfluidic gradient generators: flow-based vs diffusion-based devices (2015). *Frontiers in Bioengineering and Biotechnology*. <https://www.frontiersin.org/articles/10.3389/fbioe.2015.00039/full>
- [G4] Jeon et al. (2002). Neutrophil chemotaxis in linear and complex gradients of interleukin-8 formed in a microfabricated device. *Nature Biotechnology* 20:826–830. <https://www.gmwgroup.harvard.edu/publications/neutrophil-chemotaxis-linear-and-complex-gradients-interleukin-8-formed>

Diffusion-based separation and sensing
- [T1] Yager group tutorials on the T-sensor and H-filter, incl. measuring diffusion coefficients by fitting a numerical model (Kamholz et al. 1999, *Anal. Chem.*). <https://faculty.washington.edu/yagerp/microfluidicstutorial/tsensor/tsensor.htm>
- [T2] Weigl & Yager (1999). Microfluidic diffusion-based separation and detection. *Science*. H-filter tutorial: <https://faculty.washington.edu/yagerp/microfluidicstutorial/hfilter/hfilterhome.htm>

Oxygen and nutrients
- [O1] Consumption-driven hypoxic gradient in a cancer-cell monolayer with gas-impermeable barriers (2017). *Scientific Reports*. <https://link.springer.com/content/pdf/10.1038/s41598-017-15583-3.pdf>
- [O2] Oxygen gradients via gas-permeable PDMS membrane under cells. <https://pmc.ncbi.nlm.nih.gov/articles/PMC3188960>
- [B1] Yang et al. (2018). Analysis of factors limiting bacterial growth in PDMS mother machine devices. *Frontiers in Microbiology* 9:871. <https://pmc.ncbi.nlm.nih.gov/articles/PMC5938360>
- [B2] Wang et al. (2010). Robust growth of *Escherichia coli*. *Current Biology* (the original mother machine).

Heat
- [H1] Gui & Ren (2006). Numeric simulation of heat transfer and electrokinetic flow in an electroosmosis-based continuous-flow PCR chip. *Anal. Chem.* <https://uwaterloo.ca/waterloo-microfluidics-laboratory/sites/default/files/uploads/files/gui_ren_ac_06_0.pdf>
- [H2] Heat-loss modelling of a microfluidic PCR chip (2022). *Sensors* 22:9990. <https://mdpi.com/1424-8220/22/24/9990>

Morphogens and embryos
- [M1] Lucchetta, Lee, Fu, Patel & Ismagilov (2005). Dynamics of *Drosophila* embryonic patterning network perturbed in space and time using microfluidics. *Nature* 434:1134–1138. doi:10.1038/nature03509
- [M2] Lucchetta, Munson & Ismagilov (2006). Characterization of the local temperature in space and time around a developing *Drosophila* embryo in a microfluidic device. *Lab on a Chip*. <https://pubs.rsc.org/en/content/articlehtml/2006/lc/b516119c>
- [M3] Microfluidic embryo trap array for imaging the Dorsal gradient. <https://pmc.ncbi.nlm.nih.gov/articles/PMC4100329>
- [M4] Yu et al. (2009). Fgf8 morphogen gradient forms by a source–sink mechanism with freely diffusing molecules. *Nature* (zebrafish, FCS; not microfluidic).

Related (boundary conditions, Part 10 exercises)
- Grant et al. (2021). Predicting/simulating drug concentrations in PDMS microfluidic organ chips: drug absorption into PDMS as a partition-coefficient boundary condition. *Lab on a Chip*. <https://www.biorxiv.org/content/10.1101/2021.05.01.442248.full.pdf>

### Part 11 — Further reading (optional)

36. **`fem-outlook` — Beyond rectangles: finite elements and flow.** What
    finite volume / finite element methods do differently, meshes, why real
    chip design uses COMSOL / FEniCS, solving for the flow field itself
    (Stokes / Navier–Stokes) in outline.
37. **`iterative-solvers` — Iterative linear solvers in depth.** Jacobi,
    Gauss–Seidel, conjugate gradient, GMRES, preconditioning — for 2D/3D grids
    too large for a direct solve.

### Capstone **[?]**

A single end-to-end project after Part 9: take a time-course dataset, write an
ODE model, simulate it (Part 8), fit it (Part 9), check conditioning and
identifiability, and report uncertainty. Could double as the examination if
there is no written exam. Alternative after Part 10: design a microfluidic
gradient generator in the style of [G2] (simulate the concentration field,
then optimise the inlet concentrations so cells see a target gradient), or
recover a diffusion coefficient from a synthetic T-sensor profile [T1].

## Dependencies between parts

```
Part 1 ─┬─> Part 2 ──> (unit 7: Newton for systems) ──┬─> Part 8 (stiffness)
        │                                              └─> Part 9
        ├─> Part 3 ──> Part 4 ──> Part 5 ──> Part 8 (phase plane)
        └─> Part 6 ──> Part 7 ──> Part 8 ──> Part 10 (also needs Part 3:
                                                 sparse/tridiagonal solves)
```

Parts 2 and 3 could swap. Part 6 (interpolation) could move after Part 7 if
we want to reach ODEs sooner. **[?]**

## Open questions

- **[?]** Examination: written exam, hand-ins, capstone, or none?
- **[?]** Overlap with quantbio-3MB00A (eigenvalues, Markov chains, PCA,
  Fourier, least squares): link to those decks as refreshers, or re-teach
  from the numerical angle? Proposed: link for the concepts, teach only the
  algorithms here.
- **[?]** Part 10 can be done right after Part 8 (before optimisation) if
  the microfluidics material should come earlier. Units 31–33 are the minimum;
  34–35 could become optional.
- **[?]** Automatic differentiation: one slide (current plan) or its own unit?
