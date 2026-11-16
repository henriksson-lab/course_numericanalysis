# Cross-references from `identifiability` (unit 29)

- identifiability/lec.qmd:24 (About this lecture) -> nonlinear-least-squares (the PK fit whose parameters we now question)
- identifiability/lec.qmd:29 (About this lecture) -> svd (small singular values = badly determined combinations)
- identifiability/lec.qmd:34 (notes, "closes the loop of Part 9") -> nonlinear-least-squares, gradient-methods (fitting nonlinear models)
- identifiability/lec.qmd:67,70 (Recap: the model and the fit, "data of the previous unit") -> nonlinear-least-squares (one-compartment model, simulated data, nls/Gauss–Newton)
- identifiability/lec.qmd:116,122 (Structural (2): the flip-flop, "Recap") -> nonlinear-least-squares (flip-flop slide, multistart, algebra of V' = V ke/ka)
- identifiability/lec.qmd:135 (Near the optimum the SSE is a bowl, "Gauss–Newton approximation from the previous unit") -> nonlinear-least-squares (Gauss–Newton)
- identifiability/lec.qmd:139 (notes, "parabola from Optimisation in one dimension") -> optimisation-1d (quadratic model near a minimum)
- identifiability/lec.qmd:135-139 (Hessian) -> gradient-methods (Hessian introduced there)
- identifiability/lec.qmd:162 (notes, Standard errors from the Jacobian: numDeriv uses Richardson extrapolation) -> finite-differences (Numerical differentiation)
- identifiability/lec.qmd:221,225 (The SVD of J) -> svd (noise along v_i amplified by 1/sigma_i; CD4 - CD8 example)
- identifiability/lec.qmd:240,244 (Scaling and conditioning) -> pivoting-conditioning (lose log10 kappa digits); least-squares-qr (normal equations square kappa)
- identifiability/lec.qmd:305 (notes, Profile likelihood: uniroot) -> bisection-fixedpoint (root bracketing)
- identifiability/lec.qmd:409 (Fix 3: prior = ridge) -> svd (ridge / Tikhonov regularisation)
- identifiability/lec.qmd:413 (notes, nls.lm "from the previous unit") -> nonlinear-least-squares (Levenberg–Marquardt, minpack.lm)
- identifiability/lec.qmd:439 (Summary, Next) -> constrained-lp
- identifiability/exercise.qmd:22 (Goal) -> nonlinear-least-squares (running example and data)
- identifiability/exercise.qmd:1605 (Where this leads) -> constrained-lp (bounds ka > ke, F <= 1 as constraints)
- identifiability/exercise.qmd:1606 (Where this leads) -> capstone (contentplan "Capstone [?]", not a unit yet)
- identifiability/exercise.qmd:1607 (Where this leads) -> nonlinear-least-squares (fast, robust fitting for bootstrap/MCMC)
- identifiability/exercise.qmd: Hint 2 / B6 (eigenvalues of 2x2) -> power-iteration (eigenvalues; implicit, no explicit mention)
- identifiability/quiz.qmd:78 (Q6 hint) -> svd
- identifiability/quiz.qmd:90 (Q7 hint) -> pivoting-conditioning; least-squares-qr (normal equations)

## Running examples used here (for reconciliation in pass 2)

- **PK one-compartment oral model and data from `nonlinear-least-squares`** (checked against
  nonlinear-least-squares/make-nonlinear-least-squares-figs.R on 2026-11-16): D = 100 mg, true ka = 1 /h,
  ke = 0.2 /h, V = 20 L, noise SD 0.1 mg/L; t = (0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8, 10, 12) h,
  y = (1.02, 1.88, 2.34, 2.98, 3.27, 3.26, 3.17, 2.77, 1.92, 1.23, 1.00, 0.61) mg/L.
  Fit here (nls from (0.5, 0.1, 10)): ka 1.0150721, ke 0.1906157, V 20.2162113, SSE 0.06646744;
  flip-flop (0.1906156, 1.0150718, 3.7963082). Consistent with the nls unit's "ka = 1.02, ke = 0.19, V = 20.2".
  **Naming differs:** the nls unit calls the model `bateman(t, ka, ke, V)` with global `D` and times `tt`;
  this unit uses `pk(t, ka, ke, V, D = 100)` and `t`. Harmonise in pass 2 if wanted.

## Running examples defined here

- **Sparse design** = the same data with samples from 3 h only (`late <- t >= 3`, 6 samples):
  ka 1.1220 ± 0.4055, ke 0.1882 ± 0.0153, V 20.164 ± 0.767; profile of ka: lower bound 0.68, flat at
  ΔSSE/s² = 1.07 for large ka (no upper bound).
- **Two-compartment oral PK model** (closed-form tri-exponential, function `pk2`): ka = 1, k10 = 0.2,
  k12 = 0.3, k21 = 0.15 /h, V = 20 L; C(1), C(6), C(24) = 2.4087, 0.9606, 0.2882 mg/L; eigenvalues of
  J^T J (log-params, 12 sample times) span 4.5 decades (sloppy).
- **Enzyme product example** (exercise B1): v = kcat E0 S/(Km + S), kcat 10 /s, E0 0.1 µM, Km 2 mM.
