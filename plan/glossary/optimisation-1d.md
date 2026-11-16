**Objective function** — The function to be minimised (cost, loss), e.g. minus a log-likelihood or a sum of squared errors.
**Negative log-likelihood** — Minus the log-likelihood; minimising it gives the maximum-likelihood estimate (the convention of all optimisers).
**Stationary point** — A point where f'(x) = 0: a minimum, a maximum or a flat inflection point.
**Local minimum** — A point lower than all nearby points; not necessarily the lowest point overall.
**Global minimum** — The lowest value of f on the whole interval (or domain) considered.
**Unimodal function** — A function with a single minimum on an interval: decreasing before it, increasing after it.
**Bracket (for a minimum)** — Three points a < c < b with f(c) < f(a) and f(c) < f(b); a local minimum lies in (a, b).
**Golden-section search** — Bracketing minimiser that compares two interior points and shrinks the bracket by 0.618 per function evaluation (linear convergence).
**Golden-section ratio** — g = (√5 − 1)/2 ≈ 0.618, the factor that lets golden-section search reuse one point per step.
**Parabolic interpolation** — Fit a parabola through three points of f and jump to its vertex; superlinear (order ≈ 1.32) but unsafe far from the minimum.
**Brent's method (minimisation)** — Golden-section search combined with parabolic steps; R `optimize`, Python `scipy.optimize.minimize_scalar`, Rust argmin `BrentOpt`.
**Newton's method for minimisation** — x_{k+1} = x_k − f'(x_k)/f''(x_k): jump to the vertex of the local quadratic model; Newton's method on f'.
**Local quadratic model** — Second-order Taylor polynomial f(x_k) + f'(x_k)(x − x_k) + ½ f''(x_k)(x − x_k)^2.
**Convex function** — A function with f'' ≥ 0 everywhere on an interval; every chord lies above the graph and every local minimum is global.
**Curvature (at a minimum)** — f''(x*); a large curvature means a sharply defined minimum, for a likelihood a small standard error.
**Flat bottom (√ε limit)** — Near a minimum f changes only quadratically, so methods comparing f values locate x* only to about √ε relative.
**Profile (over linear parameters)** — Objective of the remaining nonlinear parameter after the linear ones are fitted exactly for each value, e.g. RSS(T) for a sine wave.
**Multistart** — Running a local optimiser from many starting points (or intervals) and keeping the best result, to find the global minimum.
**Grid search** — Evaluating the objective on a grid of parameter values; used to find the right basin before a local refinement.
**Half-life** — Time for a quantity decaying at rate k to halve: ln 2 / k.

## Install

No packages beyond those on intro/install.qmd (R base/stats; Python numpy, scipy, matplotlib; Rust argmin and nalgebra already in Cargo.toml).
