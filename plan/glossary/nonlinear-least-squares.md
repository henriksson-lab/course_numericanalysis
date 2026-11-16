**Nonlinear least squares** — Minimising the sum of squared residuals when the model depends nonlinearly on its parameters; solved iteratively.
**Residual sum of squares (RSS)** — $\sum_i (y_i - m_i(\mathbf p))^2$, the quantity a least-squares fit minimises (`deviance()` for `nls`).
**Sensitivity** — The change in a model prediction per unit change of a parameter; one entry of the Jacobian.
**Gauss–Newton method** — Linearise the model at the current parameters, solve the linear least-squares problem for the step (by QR), repeat.
**Levenberg–Marquardt method** — Gauss–Newton with a damping term $\lambda\|D\boldsymbol\delta\|^2$; small $\lambda$ gives Gauss–Newton, large $\lambda$ a short gradient-descent step.
**Damping parameter** — The $\lambda$ in Levenberg–Marquardt; increased after a rejected step, decreased after an accepted one.
**Marquardt scaling** — Damping with the column lengths of $J$ on the diagonal, so that the step does not depend on the units of the parameters.
**Step halving** — Shortening a Gauss–Newton step until the RSS decreases; what `nls` does.
**Singular gradient** — `nls`'s error when the Jacobian is numerically rank deficient at some trial point.
**Log-parameterisation** — Fitting $\log p$ instead of $p$, so that positive parameters stay positive and steps are relative.
**Bateman function** — The closed-form plasma concentration of the one-compartment model with first-order absorption.
**Absorption rate constant ($k_a$)** — First-order rate at which an oral drug moves from the gut into plasma.
**Elimination rate constant ($k_e$)** — First-order rate at which a drug is removed from plasma; half-life $\ln 2/k_e$.
**Volume of distribution ($V$)** — The volume that relates the amount of drug in the body to its plasma concentration.
**Flip-flop kinetics** — The symmetry $k_a \leftrightarrow k_e$ (with $V \to Vk_e/k_a$) under which oral plasma data are identical.
**Multistart** — Running a local optimiser from many starting points and comparing the results.
**Inner tolerance** — The `rtol`/`atol` of the ODE solver called inside the residual function of a fit.
**Forward sensitivity equations** — Extra ODEs whose solution is the derivative of the state with respect to a parameter; an exact alternative to finite-difference Jacobians.

## Install

- R: **minpack.lm** (`install.packages("minpack.lm")`) for `nlsLM` and `nls.lm`; it is not yet listed in `intro/install.qmd`. deSolve is also used (C5, C6).
- Python and Rust: nothing new (scipy; nalgebra, ode_solvers).
