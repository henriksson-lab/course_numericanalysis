# Cross-references from `quadrature` (Numerical integration)

- quadrature/lec.qmd:24 (About this lecture) -> intro (integrate(dnorm, 0, 20000) returning 0, exercise C5 / slide "Silent failures in two lines")
- quadrature/lec.qmd:119 (notes, "Trapezoid on a formula: halve h") -> convergence-cost (Taylor remainder, error ≈ C h^p)
- quadrature/lec.qmd:161 (notes, "Reading the order off a log–log plot") -> convergence-cost (order from the slope of a log–log plot)
- quadrature/lec.qmd:175-196 (slide "Recap: Richardson extrapolation" + notes) -> finite-differences (Richardson extrapolation of central differences, h and h/2)
- quadrature/lec.qmd:219 (notes, Romberg) -> finite-differences (difference quotients of noisy data amplify noise)
- quadrature/lec.qmd:312 (notes, "Why integrate(dnorm, 0, 20000) is 0") -> intro (the same example), floating-point (underflow of exp below the smallest double)
- quadrature/lec.qmd:335 (notes, "Three fixes") -> intro (Rust solution of C5: quadrature crate misses a peak at 5000)
- quadrature/lec.qmd:379 (slide "Biology: normalising a probability density") -> newton (zero-truncated Poisson UMI running example, MLE 1.2184695, SE 0.1325)
- quadrature/lec.qmd:391 (same slide) -> floating-point (scale by the maximum / log scale, log-sum-exp idea)
- quadrature/lec.qmd:395 (notes) -> newton (standard error used to give the scale)
- quadrature/lec.qmd:404 (slide "chromatography and mass-spec peak areas") -> nonlinear-least-squares (fitting overlapping peaks)
- quadrature/lec.qmd:450 (notes, "Data or model? AUC to infinity") -> least-squares-qr (straight-line fit of log concentration), nonlinear-least-squares (Part 9, fitting the PK model)
- quadrature/lec.qmd:476 (Summary, "Next") -> monte-carlo (integrals in many dimensions)
- quadrature/exercise.qmd:37 (Background, Richardson) -> finite-differences (Richardson extrapolation)
- quadrature/exercise.qmd:185 (C8) -> newton (zero-truncated Poisson UMI data, MLE and SE)
- quadrature/exercise.qmd:241 (Solution B5(a)) -> floating-point (smallest double, underflow)
- quadrature/exercise.qmd: Solution C7(c) -> nonlinear-least-squares (Part 9: model parameters are estimated)
- quadrature/exercise.qmd: Solution C8 -> floating-point (log-sum-exp)
- quadrature/exercise.qmd:1285 (Where this leads) -> monte-carlo (many-dimensional integrals, 1/sqrt(N))
- quadrature/exercise.qmd:1286 (Where this leads) -> euler, runge-kutta (Euler = left Riemann sum, Heun = trapezoid, adaptive step size)
- quadrature/exercise.qmd:1287 (Where this leads) -> nonlinear-least-squares (fitting the PK model and overlapping peaks)
- quadrature/quiz.qmd:70 (Q5) -> finite-differences (Richardson for difference quotients)

## Running examples defined here

- **Oral one-compartment pharmacokinetic model** (lecture, exercise, quiz, figures):
  dose 500 mg, F = 1, V = 50 L, ka = 1.5 /h, ke = 0.2 /h (t½ = 3.47 h);
  `model <- function(t) 500 / 50 * 1.5 / 1.3 * (exp(-0.2 * t) - exp(-1.5 * t))` (mg/L).
  Cmax = 7.3346 mg/L at tmax = 1.5499 h. AUC 0–6 h = 32.6243601574044,
  AUC 0–24 h = 49.5252069010181, AUC 0–inf = Dose/(V ke) = 50 mg·h/L.
  Blood samples (model rounded to 2 decimals): `t_h <- c(0, 0.5, 1, 2, 3, 4, 6, 8, 12, 24)`,
  `conc <- c(0, 4.99, 6.87, 7.16, 6.20, 5.16, 3.47, 2.33, 1.05, 0.09)`.
  Linear trapezoid AUC 0–24 = 51.6175; linear-up/log-down 48.91575; terminal slope from
  last 3 samples 0.2036783; AUC 0–inf 52.059 (linear) / 49.358 (LUL).
  Natural reuse in *Nonlinear least squares and ODE fitting* (fit ka, ke, V to the samples)
  and *Identifiability and uncertainty of fitted parameters* (ka/ke flip-flop).
- **Chromatography peak**: Gaussian, sigma = 3 s, height 1000, area 1000·3·sqrt(2π) = 7519.885;
  centred at 300 s (sampling study) or 340 s in a 0–600 s / 0–3600 s run (missed-peak failures).
