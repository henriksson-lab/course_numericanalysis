**Piecewise linear interpolation** — Joining neighbouring data points by straight lines; error at most $\tfrac{h^2}{8}\max|f''|$, never overshoots.
**Spline** — A function made of polynomial pieces joined smoothly at the data points (knots).
**Knot** — A point where two pieces of a spline meet; for interpolating and smoothing splines, the data points.
**Cubic spline** — A piecewise cubic that interpolates the data with continuous value, slope and second derivative.
**Boundary (end) condition** — One of the two extra equations a cubic spline needs, imposed at the ends of the data.
**Natural spline** — A cubic spline with zero second derivative at both ends; only $O(h^2)$ accurate near the ends.
**Clamped spline** — A cubic spline whose slopes at the two ends are prescribed.
**Not-a-knot condition** — End condition making the third derivative continuous at the second and second-to-last knot; scipy's default.
**fmm end condition** — R's default for `spline`/`splinefun`: end conditions taken from cubics through the four points at each end (Forsythe, Malcolm, Moler).
**Overshoot** — An interpolant swinging beyond the range of neighbouring data, e.g. negative expression after a steep drop.
**Monotone (shape-preserving) interpolation** — A piecewise cubic (PCHIP, Fritsch–Carlson) whose slopes are chosen so that it is monotone wherever the data are.
**Smoothing spline** — The curve minimising misfit plus $\lambda\int f''^2$; a natural cubic spline with knots at the data that does not pass through them.
**Smoothing parameter ($\lambda$)** — The weight of the roughness penalty; its size depends on the units of $t$ and $y$.
**Roughness penalty** — $\int f''(t)^2\,dt$, the total squared curvature of a curve.
**Bias–noise (bias–variance) trade-off** — Too little smoothing lets noise through; too much smoothing distorts the true curve.
**Hat matrix** — The matrix $A$ of a linear smoother, mapping the data to the fitted values, $\hat{\mathbf y} = A\mathbf y$.
**Equivalent degrees of freedom (df)** — The trace of the hat matrix; from 2 (straight line) to $n$ (interpolation); the unit-free measure of smoothing.
**Leave-one-out cross-validation (CV)** — Choosing $\lambda$ by how well each point is predicted by a fit without it; $\frac1n\sum((y_i - \hat y_i)/(1 - A_{ii}))^2$.
**Generalised cross-validation (GCV)** — $n\,\text{RSS}/(n - \text{df})^2$, CV with the leverages $A_{ii}$ replaced by their average; `smooth.spline`'s default criterion.
**Reinsch form** — Writing the smoothing spline as the banded system $(R + \lambda Q^TQ)\boldsymbol\gamma = Q^T\mathbf y$ for the second derivatives.
**Pseudotime** — An ordering of single cells along a developmental trajectory, inferred from their expression profiles.

## Install

None beyond `intro/install.qmd` (uses base R `stats`, `pracma`; scipy >= 1.10 for `make_smoothing_spline`).
