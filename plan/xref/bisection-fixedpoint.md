# Cross-references: bisection-fixedpoint

- bisection-fixedpoint/lec.qmd:33 (notes, About this lecture) -> convergence-cost (linear convergence vocabulary) [ok]
- bisection-fixedpoint/lec.qmd:33 (notes, About this lecture) -> newton (Newton, secant, uniroot/brentq come next) [ok]
- bisection-fixedpoint/lec.qmd:49 (Part 2: solving f(x)=0) -> convergence-cost (does it converge, how fast) [ok, linked]
- bisection-fixedpoint/lec.qmd:52 (notes, Part 2 slide) -> newton, newton-systems (rest of Part 2; toggle switch) [ok]
- bisection-fixedpoint/lec.qmd:172 (notes, How fast?) -> convergence-cost (linear convergence e_{k+1} = c e_k) [ok]
- bisection-fixedpoint/lec.qmd:183 (When to stop) -> floating-point (gap between doubles, eps floor) [ok, linked]
- bisection-fixedpoint/lec.qmd:209 (notes, When bisection fails) -> error-conditioning (double root is ill-conditioned) [ok]
- bisection-fixedpoint/lec.qmd:273 (notes, Why it converges) -> intro (Taylor recap), convergence-cost (linear convergence); newton (g'(x*)=0 gives faster convergence) [ok]
- bisection-fixedpoint/lec.qmd:358 (notes, Error against iteration) -> convergence-cost (rounding floor on error plots) [ok]
- bisection-fixedpoint/lec.qmd:426 (notes, Chaos slide) -> floating-point (rounding errors amplified) [ok]
- bisection-fixedpoint/lec.qmd:460 (Library routines) -> newton (how Brent's method works) [ok, linked]
- bisection-fixedpoint/lec.qmd:489 (Summary, Next) -> newton (quadratic convergence, Brent) [ok, linked]
- bisection-fixedpoint/exercise.qmd:107 (Background, library routines) -> newton (Brent's method explained) [ok, linked]
- bisection-fixedpoint/exercise.qmd:547 (Solution C2) -> optimisation-1d (maximum of surplus production = maximum sustainable catch) [ok, linked]
- bisection-fixedpoint/exercise.qmd:1094 (Where this leads) -> newton (Newton as fixed-point iteration with g'(x*)=0; Brent) [ok, linked]
- bisection-fixedpoint/exercise.qmd:1095 (Where this leads) -> newton-systems (toggle switch steady states) [ok, linked]
- bisection-fixedpoint/exercise.qmd:1096 (Where this leads) -> phase-plane (continuous-time stability: Jacobian eigenvalues) [ok, linked]
- bisection-fixedpoint/exercise.qmd:1097 (Where this leads) -> optimisation-1d (golden-section search for the maximum sustainable catch) [fixed: optimisation-1d does not use the fish/MSY example; reworded to "a maximum golden-section search would find (that unit uses its own examples)"]

## Running examples defined here

- **Harvested Ricker fish stock** (thousands of fish): N_{t+1} = g(N) = N exp(r (1 - N/K)) - H with r = 0.5, K = 100, H = 10.
  Steady states (roots of f(N) = g(N) - N): 20.4819903614 (unstable, g' = 1.336, collapse threshold)
  and 74.9524424994 (stable, g' = 0.709). Brackets [0, 50] and [50, 100]. Maximum sustainable catch
  (maximum of surplus N(exp(r(1-N/K)) - 1)) = 14.26 at N = 46.75; above it no steady state.
  Rearrangements g2(N) = K(1 - log(1 + H/N)/r) (finds upper root, rate 0.314) and
  g3(N) = H/(exp(r(1-N/K)) - 1) (finds lower root, rate 0.312). Useful for *Newton and secant methods*
  (compare iteration counts: bisection 33 steps to 1e-8, uniroot 9 to 1e-12) and *Optimisation in one dimension* (MSY).
- **Ricker without harvest**, K = 100, N0 = 20: g'(K) = 1 - r; r = 0.5, 1.5 (converge), 2.3 (2-cycle 159.22/40.78),
  2.6 (4-cycle), 3.0 (chaos). Period doubling at r = 2, 2.53, chaos beyond 2.69, 3-cycle window near 3.12.
- **Beverton–Holt**, R = 2, K = 100, N0 = 20 (slides: 33.33, 50, 66.67, 80, ...); exercise uses R = 3, N0 = 50; g'(K) = 1/R.
