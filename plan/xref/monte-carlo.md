# Cross-references from `monte-carlo` (unit 20)

- monte-carlo/lec.qmd:24 (About this lecture) -> quadrature (grid / tensor-product quadrature rules) [ok]
- monte-carlo/lec.qmd:34 (notes, About this lecture, "quadrature rules from the previous unit") -> quadrature (high-order quadrature rules) [ok]
- monte-carlo/lec.qmd:104 (notes, A second problem: an integral in 10 dimensions, "tensor-product rules of the previous unit") -> quadrature (cost n^d of tensor-product rules) [ok]
- monte-carlo/lec.qmd:204 (notes, Inverse-transform sampling) -> gillespie (exponential waiting times by inverse transform) [ok]
- monte-carlo/lec.qmd:258 (notes, Monte Carlo integration) -> quadrature (AUC of the OD curve with the trapezoid rule) [fixed: quadrature computed a PK AUC, not the OD AUC; notes now say 'the same kind of AUC ... (there for a drug concentration curve)']
- monte-carlo/lec.qmd:258 (notes, Monte Carlo integration) -> convergence-cost (shared logistic OD running example K = 1.2, N0 = 0.05, r = 0.8) [ok]
- monte-carlo/lec.qmd:291 (notes, The error falls like 1/sqrt(N)) -> quadrature (trapezoid rule, order 2; integrate/quad) [ok]
- monte-carlo/lec.qmd:312 (The curse of dimensionality, table "trapezoid-like rule of order 2") -> quadrature (order of composite rules) [ok]
- monte-carlo/lec.qmd:320 (notes, The curse of dimensionality) -> quadrature (tensor-product rules, Simpson order 4, Gauss rules) [ok]
- monte-carlo/lec.qmd:344 (notes, The marginal likelihood in 10 dimensions, "floating-point unit") -> floating-point (log-likelihoods avoid underflow) [ok]
- monte-carlo/lec.qmd:441 (Summary, Next) -> euler (start of Part 8) [ok]
- monte-carlo/lec.qmd:441 (Summary) -> gillespie (random numbers drive Gillespie's algorithm) [ok]
- monte-carlo/lec.qmd:444 (notes, Summary) -> gillespie (inverse-transform sampling of exponential waiting times) [ok]
- monte-carlo/exercise.qmd:159 (C4, "the running example") -> convergence-cost (logistic OD curve K = 1.2, N0 = 0.05, r = 0.8) [ok]
- monte-carlo/exercise.qmd:Hint 4 (log-scale likelihood to avoid underflow) -> floating-point (log-likelihoods) [ok]
- monte-carlo/exercise.qmd:1144 (Where this leads) -> gillespie (exponential waiting times, discrete inverse transform for the next reaction) [ok: gillespie uses inverse transform for exponential waiting times and the discrete inverse transform (findInterval on cumsum), and SE of averages]
- monte-carlo/exercise.qmd:1145 (Where this leads) -> identifiability (bootstrap as a Monte Carlo average) [fixed: identifiability teaches SEs and profiles, only mentions the bootstrap in Where this leads; wording changed]

## Running examples defined here

- **Motif hits in promoters** (lec slides 3–5, exercise B6, C7, quiz Q8): region G = 1e6 bp, 50 promoters of
  1000 bp starting at 10001 + 20000 (j − 1), j = 1..50; motif length L = 10; 200 hits, 18 observed in promoters.
  p_hit = 50·1009/999991 = 0.05045, E[K] = 10.09, exact P(K ≥ 18) = 0.01314; R set.seed(1), 10 000 placements
  with sample.int → 122 (0.0122). Exact P(K ≥ 25) = 3.0e-5.
- **Toy marginal likelihood** (lec, exercise C6): d parameters, prior N(0,1), one measurement y_i = 1 with noise
  sd 1; Z = dnorm(1, 0, sqrt(2))^d (0.2196956 per dimension; 2.619476e-07 for d = 10). Grid at 3 prior quantiles per
  axis in d = 10 is 78% wrong. Relative variance of L under the prior 1.364^d − 1.
- **Rare motif score** P(X > 5) = 2.866516e-07 for X ~ N(0,1), importance sampling with N(5,1) (weights exp(12.5 − 5x)).
- **Burst waiting times** exponential with rate 0.5 per hour (inverse transform) — could be reused in *Stochastic simulation*.
- OD-curve AUC on [0, 10] h for the shared logistic example: exact (K/r) log((e^{8} + 23)/24) = 7.244448.

## Pass 2

- All italic unit titles turned into links (first mention per slide / page / Where this leads).
