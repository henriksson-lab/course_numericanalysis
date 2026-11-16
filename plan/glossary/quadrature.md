**Quadrature** — Numerical computation of an integral as a weighted sum of function values, ∑ w_i f(x_i).
**Area under the curve (AUC)** — Integral of drug concentration over time (mg·h/L); the standard measure of total drug exposure.
**Riemann sum** — Approximation of an integral by rectangles of width h; the left-endpoint version has order 1.
**Midpoint rule** — Rectangles with height f at the centre of each interval; order 2, never evaluates f at the endpoints.
**Trapezoid rule** — Straight line between neighbouring points, area h(f_{i-1} + f_i)/2 per interval; order 2.
**Simpson's rule** — Parabola through three equally spaced points, weights h/3 (1, 4, 1); order 4, exact for cubics.
**Composite rule** — A simple rule applied on each of n subintervals and summed.
**Order of a quadrature rule** — The power p in error ≈ C h^p; minus the slope of error against n on a log–log plot.
**Euler–Maclaurin formula** — Expansion of the trapezoid error in even powers h², h⁴, … with coefficients from derivatives at the endpoints.
**Romberg integration** — Repeated Richardson extrapolation of trapezoid sums with h, h/2, h/4, …, removing h², h⁴, … in turn.
**Adaptive quadrature** — Estimate the error on each piece by comparing two rules and subdivide only pieces that fail the tolerance.
**Adaptive Simpson** — Adaptive quadrature comparing Simpson on one panel and on two half-panels; |S2 − S1|/15 estimates the error.
**Gaussian (Gauss–Legendre) quadrature** — Rule with optimally chosen points and weights; n points integrate polynomials of degree 2n − 1 exactly.
**Gauss–Kronrod rule** — A Gauss rule extended with extra points so the difference of the two gives an error estimate; QUADPACK uses 10 + 11 = 21 points.
**QUADPACK** — Fortran library of adaptive quadrature routines (QAGS, QAGI) behind R's integrate and SciPy's quad.
**Double-exponential (tanh–sinh) rule** — Quadrature after a change of variables that packs points towards both ends; robust to endpoint singularities (Rust quadrature crate).
**Error estimate** — The routine's own guess of its error, computed from the same points as the value; can be 0 when every point misses a feature.
**Absolute / relative tolerance** — Stopping thresholds of a quadrature routine; an integral smaller than the absolute tolerance is accepted after one rough estimate.
**Integrable singularity** — A point where the integrand is infinite but the area is finite, e.g. 1/√x at 0.
**Substitution (change of variables)** — Rewrite the integral in a new variable (x = u², x = t/(1 − t), x = sin²θ) to remove a singularity or an infinite interval.
**Normalising constant** — The integral Z of an unnormalised density, so that density/Z integrates to 1; the marginal likelihood in Bayesian work.
**Jeffreys prior** — For a proportion, the Beta(½, ½) density ∝ 1/√(x(1 − x)); normalising constant π.
**Linear-up / log-down rule** — AUC rule using the trapezoid while concentrations rise and the area of an exponential, Δt (C1 − C2)/ln(C1/C2), while they fall.
**Terminal slope (λz)** — Elimination rate from a straight-line fit of log concentration over the last samples; used to extrapolate AUC to infinity as C_last/λz.
**Non-compartmental analysis (NCA)** — Pharmacokinetic analysis directly from sampled concentrations (AUC by quadrature), without fitting a model.
**Bioequivalence** — Two formulations are judged equivalent when the ratio of their AUCs (and Cmax) lies within 80–125 %.
**Peak area** — Integral of a chromatography or mass-spectrometry signal over a peak; the quantity used for quantification, more robust than peak height.

## Install

None beyond `intro/install.qmd` (R: base/stats; pracma only named on a slide for Gauss–Legendre points, not needed to run anything; Python: numpy, scipy; Rust: quadrature, statrs, already in `rust/Cargo.toml`).
