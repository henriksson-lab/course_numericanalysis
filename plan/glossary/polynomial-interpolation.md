**Interpolation** — Constructing a function that passes exactly through given data points, used to estimate values between them.
**Interpolating polynomial** — The unique polynomial of degree at most $n$ through $n+1$ points with distinct $x$-values.
**Nodes** — The points $t_0, \dots, t_n$ (e.g. sampling times) at which the data are given.
**Vandermonde matrix** — The matrix with entries $t_i^{\,j}$; solving $V\mathbf c = \mathbf y$ gives the monomial coefficients of the interpolant, but it is badly conditioned.
**Lagrange form** — The interpolant written as $\sum_i y_i L_i(t)$, a weighted sum of the readings.
**Lagrange basis polynomial** — $L_i(t) = \prod_{j\ne i}(t - t_j)/(t_i - t_j)$: 1 at node $t_i$, 0 at the other nodes.
**Newton form** — The interpolant written as $a_0 + a_1(t-t_0) + a_2(t-t_0)(t-t_1) + \dots$; adding a node adds one term.
**Divided difference** — The recursively defined quotient $f[t_0, \dots, t_k]$; the Newton coefficients, equal to $f^{(k)}(\xi)/k!$ for some $\xi$.
**Interpolation error formula** — $f(t) - p(t) = f^{(n+1)}(\xi)\,w(t)/(n+1)!$: a derivative of the curve times the node polynomial.
**Node polynomial** — $w(t) = \prod_i (t - t_i)$; depends only on where the data were sampled.
**Runge phenomenon** — Growing oscillations near the ends of the interval when a high-degree polynomial interpolates at equally spaced nodes.
**Chebyshev nodes** — Points $\cos((2k+1)\pi/(2n+2))$ (mapped to $[a, b]$), crowded at the ends; they minimise $\max|w(t)|$ and avoid the Runge phenomenon.
**Lebesgue constant** — $\max_t \sum_i |L_i(t)|$, the worst-case factor by which errors in the data are amplified by interpolation (89 for 13 equispaced nodes, 2.6 for 13 Chebyshev nodes).
**Barycentric formula** — $p(t) = \sum_j \frac{w_j y_j}{t - t_j} \big/ \sum_j \frac{w_j}{t - t_j}$; the stable $O(n)$-per-point way to evaluate an interpolating polynomial.
**Barycentric weights** — $w_j = 1/\prod_{k\ne j}(t_j - t_k)$; only their ratios matter.
**Piecewise linear interpolation** — Joining neighbouring readings by straight lines (R `approx`, `np.interp`); error at most $h^2\max|f''|/8$.
**Extrapolation** — Evaluating an interpolant or fit outside the range of the data; for polynomials, unreliable.
**Interpolation vs fitting** — Interpolation passes through every point (as many parameters as points); fitting (least squares) uses fewer parameters and averages out noise.

## Install

No packages beyond those in `intro/install.qmd` (base R `lm`, `poly`, `approx`, `solve`, `kappa`; numpy `numpy.polynomial`, scipy `scipy.interpolate.BarycentricInterpolator`; nalgebra in Rust).
