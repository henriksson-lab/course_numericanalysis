**Identifiable parameter** — a parameter whose value is determined by the model output: different values give different predictions.
**Structural non-identifiability** — the model equations give identical output for different parameter values (e.g. F and V only as F/V); no data can fix it.
**Practical non-identifiability** — the model could determine the parameter, but the available data (too few, too noisy, badly timed) do not.
**Flip-flop** — the symmetry of the oral one-compartment model: (ka, ke, V) and (ke, ka, V·ke/ka) give the same concentration curve.
**Apparent volume (V/F)** — the identifiable combination of volume and bioavailability after an oral dose.
**Standard error (of a fitted parameter)** — square root of a diagonal element of the approximate covariance matrix s²(JᵀJ)⁻¹.
**Wald interval** — estimate ± 1.96 standard errors; relies on the sum of squares being a quadratic bowl.
**Fisher information** — JᵀJ/σ² for least squares with Gaussian noise (in general the Hessian of −log L); its inverse approximates the covariance of the estimates.
**Correlation matrix of estimates** — covariance scaled by the standard errors; |r| near 1 means two parameters can compensate each other.
**Sensitivity (relative)** — θⱼ ∂C/∂θⱼ, the change in the output per relative change in a parameter; a column of the log-parameter Jacobian.
**Profile likelihood** — the minimal sum of squares (or −log L) over all other parameters with one parameter held fixed; its rise above the minimum gives a confidence interval.
**Flat profile** — a profile that stays below the 3.84 threshold in one or both directions: the parameter is not bounded there by the data.
**Sloppy model** — a model whose JᵀJ eigenvalues (in log-parameters) are spread over many decades without a gap: a few stiff and many poorly determined parameter combinations.
**Delta method** — propagating the parameter covariance to a derived quantity g: Var g ≈ ∇gᵀ Cov ∇g.
**Optimal experimental design** — choosing sample times (or doses) to make JᵀJ large in the directions of interest, before the experiment.
**Prior as a residual** — adding (log θ − μ)/τ to the residual vector; equivalent to ridge regularisation and to the maximum a posteriori estimate with a log-normal prior.

## Install

- R: `numDeriv` (jacobian, hessian) and `minpack.lm` (nls.lm); neither is on intro/install.qmd yet (earlier units also flagged numDeriv).
- Python: nothing beyond numpy/scipy (uses `scipy.differentiate.hessian`, scipy ≥ 1.15).
