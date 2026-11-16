# Cross-references from `optimisation-1d` (Optimisation in one dimension)

- optimisation-1d/lec.qmd:34 (notes, About this lecture) -> gradient-methods (next unit, several parameters); bisection-fixedpoint, newton (method twins: bisection/golden, secant/parabolic, Brent/Brent, Newton/Newton)
- optimisation-1d/lec.qmd:37-47 (slide "Recap: an MLE as the root of the score") -> newton (zero-truncated Poisson MLE, score, Newton from the mean, uniroot)
- optimisation-1d/lec.qmd:173 (notes, Why 0.618?) -> convergence-cost (linear convergence vocabulary)
- optimisation-1d/lec.qmd:170 (Why 0.618?, last bullet) -> bisection-fixedpoint (bisection 0.30 digits per step)
- optimisation-1d/lec.qmd:229 (notes, Parabolic interpolation) -> newton (secant method as twin)
- optimisation-1d/lec.qmd:240-257 (slide "Brent's method") -> newton (uniroot / safeguarded Brent root finder)
- optimisation-1d/lec.qmd:285 (notes, optimize in practice) -> newton (uniroot default tol, <<- counter from the Newton exercise C6)
- optimisation-1d/lec.qmd:290-320 (slides "How precisely..." and "The flat bottom") -> floating-point (machine epsilon, gap between doubles), newton (score, expm1)
- optimisation-1d/lec.qmd:329 (Newton for minimisation) -> newton (Newton on the score, 6 steps from the mean)
- optimisation-1d/lec.qmd:334 (notes) -> identifiability (standard errors from curvature)
- optimisation-1d/lec.qmd:370 (notes, When Newton heads uphill) -> gradient-methods (line search, safeguards)
- optimisation-1d/lec.qmd:407 (notes, circadian period) -> least-squares-qr (lm.fit for the linear parameters), fourier-fft (periodogram)
- optimisation-1d/lec.qmd:443 (notes, Grid first) -> fourier-fft (periodogram = grid scan), gradient-methods / nonlinear-least-squares (multistart, starting values)
- optimisation-1d/lec.qmd:456 (notes, Checklist) -> constrained-lp (optimum on a boundary)
- optimisation-1d/lec.qmd:467 (Summary, Next) -> gradient-methods
- optimisation-1d/lec.qmd:470 (notes, Summary) -> gradient-methods (line search), nonlinear-least-squares (two-parameter decay fit), identifiability (curvature -> standard errors)
- optimisation-1d/exercise.qmd:44 (Background, running examples) -> newton (zero-truncated Poisson data)
- optimisation-1d/exercise.qmd:183 (C3(c)) -> newton (S'(lambda) formula)
- optimisation-1d/exercise.qmd:196 (C6) -> convergence-cost (logistic OD curve K = 1.2, N0 = 0.05, r = 0.8), finite-differences (the 61 noisy OD values)
- optimisation-1d/exercise.qmd:261 (Solution B4(e)) -> gradient-methods (line search)
- optimisation-1d/exercise.qmd:1047 (Solution C6) -> finite-differences (3-parameter MLE r = 0.7954), gradient-methods
- optimisation-1d/exercise.qmd:1053-1056 (Where this leads) -> gradient-methods, nonlinear-least-squares, identifiability, fourier-fft
- optimisation-1d/exercise.qmd: Hint 6 (counting calls) -> newton (AtomicU64 counter as in Newton exercise C6)
- optimisation-1d/quiz.qmd:142 (Q11) -> constrained-lp (optimum on a boundary)

## Running examples defined here

- **mRNA decay time course** (fraction left after transcription stop, qPCR; simulated k = 0.35/h, sd 0.03, seed 7):
  t = 0, 0.5, 1, 2, 4, 6, 8 h; y = 1.000, 0.804, 0.684, 0.484, 0.217, 0.094, 0.083.
  Model exp(-k t), SSE(k) = sum (y - exp(-k t))^2. k-hat = 0.3753873095 /h (root of SSE'), SSE = 0.002042809,
  half-life 1.846 h. SSE'' < 0 for k > 0.8021 (inflection). log-linear lm gives 0.3476 (different error model).
  Newton on SSE' from 0.2 converges in 7 steps; from 0.7 jumps to -1.66; from 1.0 runs to +infinity.
  optimize default tol: 0.3754002 in 11 evals; tol 1e-10: 14 evals; golden 1e-10: 50 evals.
  Natural reuse in nonlinear-least-squares (two parameters: initial amount and rate).
- **Circadian bioluminescence** (simulated: 2 + 0.6 cos(2 pi (t - 6)/24.5) + N(0, 0.3^2), set.seed(15), rounded to 2 dp),
  t = 0, 4, ..., 140 h (36 points); values hard-coded in exercise C5.
  RSS(T) of the best sine wave (lm.fit on 1, cos, sin): global min T = 24.63047 h (RSS 1.908081),
  local min T = 32.40162 h (RSS 10.225063); RSS of a constant 10.72143. optimize on [10, 40] -> 32.40 (first points
  21.46, 28.54), on [15, 35] -> 24.63. Dip width about 3 h (RSS < 6 for 23.0-26.25 h).
- **ZTP flat bottom**: nll(lhat) = 52.61033, nll'' = 56.9186; predicted indistinguishable band 2.03e-8;
  optimize(tol = 1e-15) error 6.86e-10 (R), 6.71e-10 (scipy bounded), 2.18e-9 (argmin BrentOpt).
- **OD growth rate with K, N0 fixed** (data of finite-differences): r-hat = 0.800270322, SSE = 0.001232394, residual sd 0.00453.
