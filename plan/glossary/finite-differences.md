**Finite difference** — an estimate of a derivative from function values at nearby points, divided by a power of the step h.
**Stencil** — the set of offsets and weights of a finite-difference formula, e.g. (1, −2, 1)/h² for the second derivative.
**Forward difference** — (f(x+h) − f(x))/h; first-order estimate of f′, error (h/2) f″.
**Backward difference** — (f(x) − f(x−h))/h; first-order estimate of f′ using only past values.
**Central difference** — (f(x+h) − f(x−h))/(2h); second-order estimate of f′, error (h²/6) f‴.
**One-sided difference** — a stencil using points on one side only, e.g. (−3f(x) + 4f(x+h) − f(x+2h))/(2h) at the start of a data set.
**Second difference** — (f(x+h) − 2f(x) + f(x−h))/h²; second-order estimate of f″, the stencil of the diffusion term.
**Method of undetermined coefficients** — finding stencil weights by matching Taylor terms, a small linear system.
**Optimal step size** — the h that balances truncation error C hᵖ against rounding error ε/hᵏ; about ε^(1/(p+k)).
**Richardson extrapolation** — combining estimates at h and h/2 to cancel the leading error term, e.g. (4D(h/2) − D(h))/3.
**Richardson table** — repeated Richardson extrapolation, each column one order higher; used by numDeriv::grad.
**Noise amplification** — differentiating measured data divides the measurement noise by the step, so the derivative is noisier than the data.
**Smoothing before differentiating** — fitting a smooth curve (local fit, smoothing spline) to noisy data and differentiating the curve.
**Savitzky–Golay filter** — a local least-squares polynomial fit in a moving window, used to smooth or differentiate sampled data.
**Gradient** — the vector of partial derivatives of a scalar function, e.g. of a log-likelihood with respect to the parameters.
**Gradient check** — comparing a hand-written or AD gradient with a finite-difference gradient to find bugs.
**Complex-step derivative** — Im f(x + ih)/h; a derivative without subtraction, accurate to rounding for code that accepts complex numbers.
**Automatic differentiation (AD)** — applying the chain rule to every elementary operation of a program to get derivatives exact to rounding.
**Forward-mode AD** — AD that carries (value, derivative) pairs through the computation; one pass per input direction.
**Reverse-mode AD** — AD that records the computation and applies the chain rule backwards; the whole gradient for a few function evaluations (backpropagation).

## Install

- R: `numDeriv` (numDeriv::grad, numDeriv::hessian). Not yet on intro/install.qmd (newton-systems also needs it).
- Python: nothing new (scipy ≥ 1.15 for `scipy.differentiate.derivative`; `scipy.interpolate.make_smoothing_spline` needs scipy ≥ 1.10).
- Rust: nothing new (nalgebra only).
