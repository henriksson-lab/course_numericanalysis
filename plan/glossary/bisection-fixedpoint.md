**Root** — a number x* with f(x*) = 0.
**Bracket** — an interval [a, b] where f(a) and f(b) have opposite signs; for continuous f it contains a root.
**Intermediate value theorem** — a continuous function that changes sign on an interval has a root in it.
**Bisection** — root finding by halving a bracket and keeping the half with the sign change; linear convergence with rate 1/2.
**Double root** — a root where f touches zero without changing sign (f and f' both zero); invisible to bracketing methods.
**Pole** — a point where f becomes infinite; it can produce a sign change without a root.
**Residual** — the value f(x) at a computed root; should be close to zero.
**Stopping criterion** — the rule that ends an iteration: bracket width or step size below an absolute plus relative tolerance, or a maximum number of iterations.
**Fixed point** — a number x* with g(x*) = x*.
**Fixed-point iteration** — repeating x_{k+1} = g(x_k) to find a fixed point of g.
**Cobweb plot** — picture of a fixed-point iteration: alternate vertical steps to y = g(x) and horizontal steps to the diagonal y = x.
**Attracting (stable) fixed point** — a fixed point with |g'(x*)| < 1; nearby iterates converge to it.
**Repelling (unstable) fixed point** — a fixed point with |g'(x*)| > 1; nearby iterates move away from it.
**Rate of linear convergence** — the factor L = |g'(x*)| by which the error shrinks each step.
**Steady state (equilibrium)** — a population size that does not change from one generation to the next, N* = g(N*).
**Beverton–Holt model** — discrete-time population model N_{t+1} = R N_t / (1 + (R-1) N_t / K); always approaches K monotonically.
**Ricker model** — discrete-time population model N_{t+1} = N_t exp(r (1 - N_t/K)); overshoots, cycles or becomes chaotic for large r.
**Maximum sustainable yield** — the largest constant harvest for which a steady state exists.
**Cycle (period-k orbit)** — a set of k values that an iteration visits in turn forever, e.g. a 2-cycle of boom and bust years.
**Bifurcation diagram** — plot of the long-run values of an iteration against a parameter.
**Chaos** — bounded, non-periodic behaviour of a deterministic iteration with sensitive dependence on the starting value.
**Sensitive dependence on initial conditions** — tiny differences in the start (or rounding errors) grow exponentially.
**Brent's method** — safeguarded root finder used by uniroot and brentq: keeps a bracket like bisection but takes faster interpolation steps.

## Install

None beyond intro/install.qmd (R: base/stats; Python: numpy, scipy, matplotlib; Rust: argmin, already in Cargo.toml).
