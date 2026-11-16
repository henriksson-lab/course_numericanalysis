# Cross-references from polynomial-interpolation

- polynomial-interpolation/lec.qmd:33 -> splines (About: practical alternative is the next lecture; notes) [ok]
- polynomial-interpolation/lec.qmd:45 -> convergence-cost (running logistic OD curve K = 1.2, N0 = 0.05, r = 0.8; notes) [fixed: slide now links Convergence and cost; K, N0, r agree]
- polynomial-interpolation/lec.qmd:67 -> least-squares-qr (polynomial fit = Vandermonde system; interpolation = as many coefficients as points) [fixed: linked]
- polynomial-interpolation/lec.qmd:70 -> least-squares-qr (Vandermonde matrix, conditioning; notes) [ok]
- polynomial-interpolation/lec.qmd:138 -> convergence-cost (first divided difference = forward-difference slope; notes) [ok]
- polynomial-interpolation/lec.qmd:197 -> convergence-cost (Taylor remainder, same shape as interpolation error) [fixed: linked]
- polynomial-interpolation/lec.qmd:222 -> convergence-cost (f''(3) = 0.0592 set the forward-difference error constant; notes) [ok: convergence-cost has N''(3) = 0.0592, constant 0.0296]
- polynomial-interpolation/lec.qmd:287 -> splines ("cubic splines of the next lecture") [fixed: linked]
- polynomial-interpolation/lec.qmd:295 -> splines (keep locality and add smoothness; notes) [ok]
- polynomial-interpolation/lec.qmd:306 -> splines (bias-variance, smoothing splines; notes) [ok]
- polynomial-interpolation/lec.qmd:330 -> least-squares-qr (lm with fewer coefficients is least squares; poly() is a QR of the Vandermonde matrix; notes) [ok]
- polynomial-interpolation/lec.qmd:348 -> least-squares-qr (Vandermonde condition number grows exponentially; notes) [ok]
- polynomial-interpolation/lec.qmd:380 -> splines (checklist: smooth noisy readings) [fixed: linked]
- polynomial-interpolation/lec.qmd:397 -> splines (Next: slide) [fixed: linked]
- polynomial-interpolation/lec.qmd:400 -> splines (growth rate as derivative of smoothed OD curve; notes) [ok]
- polynomial-interpolation/exercise.qmd:231 -> least-squares-qr (Hint 5: lose log10 kappa digits) [fixed: linked (first mention on page)]
- polynomial-interpolation/exercise.qmd:647 -> convergence-cost (C2 solution: f''(3) = 0.0592) [ok]
- polynomial-interpolation/exercise.qmd:1446 -> splines (Where this leads: local error, no Runge, smoothing splines, growth rate) [fixed: linked]
- polynomial-interpolation/exercise.qmd:1447 -> fourier-fft (Where this leads: trigonometric interpolation, circadian data) [fixed: linked]
- polynomial-interpolation/exercise.qmd:1448 -> finite-differences, quadrature (Where this leads: differences/Simpson from interpolating polynomials; Gaussian quadrature chooses nodes) [fixed: linked; reworded: finite-differences builds stencils from Taylor series, not from interpolating polynomials; quadrature does cover Simpson-as-parabola and Gaussian quadrature]
- polynomial-interpolation/exercise.qmd:1449 -> least-squares-qr (Where this leads: polynomial fits and conditioning) [fixed: linked]
- polynomial-interpolation/quiz.qmd:Q7 -> (internal: "More readings make it worse" slide of this lecture) [ok]

## Running examples defined here

- **48-h plate-reader run** (useful for `splines`): the course logistic OD curve read every 4 h, `t48 <- seq(0, 48, by = 4)`; the degree-12 interpolant dips to OD -0.148 at 0.87 h (max error 0.252); at 13 Chebyshev times on [0, 48] max error 0.0386.
- **Noisy hourly OD readings** (useful for `splines` smoothing, `finite-differences` growth rate from data): `t12 <- 0:12`, `y12 <- c(0.0375, 0.1095, 0.1959, 0.4207, 0.6259, 0.8279, 1.0188, 1.1207, 1.1673, 1.1736, 1.2210, 1.2037, 1.1857)` (logistic curve + N(0, 0.02^2), `set.seed(1)`, rounded to 4 decimals). Degree-12 interpolant max error 0.394; least-squares degree 4 max error 0.030 (residual sd 0.023); `approx` max error 0.032.
- **Immediate-early gene time course** (useful for `splines`, `fourier-fft` contrast): fold change `g(t) = 1 + 6 * (exp(-0.25 * t) - exp(-1.5 * t))` sampled at `c(0, 0.5, 1, 2, 4, 8, 12, 24)` h; degree-7 interpolant reaches -2241 at 21.4 h, -1120 at 18 h; `approx` max error 0.283.

## Pass 2 notes

- OD curve parameters and f''(3) = 0.0592, max|f''| = 0.0739 agree with convergence-cost and splines.
- "needs ... from C3" code comments now say the functions are defined in the C3 solution on this page.
