**Newton's method** — Root-finding iteration x_{k+1} = x_k - f(x_k)/f'(x_k): solve the tangent line instead of f.
**Tangent line** — First-order Taylor polynomial f(x_k) + f'(x_k)(x - x_k); Newton's local model of f.
**Secant method** — Newton's method with f' replaced by the slope of the chord through the last two iterates; no derivative needed.
**Quadratic convergence** — Error squared each step, e_{k+1} ≈ C e_k^2; the number of correct digits roughly doubles (order 2).
**Superlinear convergence** — Faster than linear (order p > 1) but not necessarily quadratic; the secant method has p ≈ 1.618.
**Order of convergence** — The power p in e_{k+1} ≈ C e_k^p; estimated as log(e_{k+1}/e_k) / log(e_k/e_{k-1}).
**Golden ratio** — (1 + √5)/2 ≈ 1.618, the order of convergence of the secant method.
**Simple root** — A root x* with f'(x*) ≠ 0; Newton converges quadratically near it.
**Double (multiple) root** — A root where f'(x*) = 0 as well, e.g. (x - 1)^2; Newton converges only linearly and no sign change exists.
**Score function** — Derivative of the log-likelihood with respect to a parameter; the maximum-likelihood estimate is a root of it.
**Observed information** — Minus the second derivative of the log-likelihood at the MLE, -S'(θ̂); 1/√ of it is the usual standard error.
**Zero-truncated Poisson** — Poisson distribution conditioned on a count ≥ 1, for data where zeros are not observed.
**EC50 / c50** — Dose that gives 50 % of the response; for a single Hill curve it equals K, for mixtures it must be found by root finding.
**Two-site (biphasic) dose–response** — Response that is a weighted sum of two Hill curves, e.g. from two receptor subtypes.
**Cycling** — An iteration that returns to earlier values and repeats forever instead of converging (e.g. Newton on a sigmoid from a bad start).
**Safeguarded method** — Root finder that keeps a sign-change bracket and falls back to bisection when a fast step would leave it.
**Brent's method** — Safeguarded root finder combining bisection, secant and inverse quadratic interpolation; R `uniroot`, Python `scipy.optimize.brentq`.
**Inverse quadratic interpolation** — Fit x as a quadratic function of f through three points and evaluate at f = 0; the fast step in Brent's method.
**Function evaluation count** — Number of calls of f (and f'); the fair measure of a root finder's cost.

## Install

None beyond `intro/install.qmd` (R: base/stats and numDeriv only used in a teacher check; Python: numpy, scipy; Rust: argmin, already in `rust/Cargo.toml`).
