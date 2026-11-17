**Brownian motion (Wiener process)** — $W(t)$ with $W(0) = 0$ and independent increments $W(t + h) - W(t) \sim N(0, h)$; continuous but nowhere differentiable.
**Stochastic differential equation (SDE)** — $dX = f(X)\,dt + g(X)\,dW$: an ODE (drift $f$) plus noise of strength $g$ driven by Brownian motion.
**Drift** — the deterministic part $f$ of an SDE; the right-hand side of the underlying ODE.
**Diffusion coefficient (of an SDE)** — the noise strength $g$; in the Fokker–Planck equation it enters as $D = g^2/2$.
**Additive noise** — noise strength $g$ that does not depend on $X$.
**Multiplicative noise** — noise strength $g(X)$ that depends on $X$, e.g. $\sigma N$ for a population.
**Euler–Maruyama method** — $X_{n+1} = X_n + f(X_n)h + g(X_n)\sqrt h\,Z_n$; strong order ½, weak order 1 for multiplicative noise.
**Milstein method** — Euler–Maruyama plus $\tfrac12 g g'((\Delta W)^2 - h)$; strong order 1 for scalar noise.
**Strong error** — mean distance between numerical and exact paths driven by the same Brownian path.
**Weak error** — error in an average (mean, variance, probability) of the numerical solution.
**Mean squared displacement (MSD)** — average of $|\mathbf r(t) - \mathbf r(0)|^2$; equals $2dDt$ for free diffusion in $d$ dimensions.
**Stokes–Einstein relation** — $D = k_BT/(6\pi\eta a)$ for a sphere of radius $a$ in a liquid of viscosity $\eta$.
**Ornstein–Uhlenbeck process** — $dX = -\theta X\,dt + \sigma\,dW$; mean-reverting, stationary variance $\sigma^2/(2\theta)$.
**Geometric Brownian motion** — $dN = rN\,dt + \sigma N\,dW$; mean grows at $r$, median at $r - \sigma^2/2$.
**Itô correction** — the term $-\sigma^2/2$ in the growth rate of $\log N$, from $(\Delta W)^2 \approx h$.
**Itô vs Stratonovich** — two interpretations of the noise term (evaluated at the start or the midpoint of a step); Euler–Maruyama gives Itô.
**Mean-square stability** — $\mathbf E\,X_n^2 \to 0$ for a decaying linear SDE; Euler–Maruyama may need a smaller step for it than Euler for the ODE.
**Drift-implicit Euler–Maruyama** — Euler–Maruyama with the drift evaluated at the new point (backward Euler for the drift).
**Fokker–Planck equation** — PDE for the probability density of an SDE's solution: advection with velocity $f$, diffusion with $D = g^2/2$.
**Common random numbers** — using the same random increments in two simulations to compare them with little Monte Carlo noise.

## Install

None (base R, numpy, scipy; Rust rand, rand_distr, statrs).
