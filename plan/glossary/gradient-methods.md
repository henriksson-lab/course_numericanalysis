**Gradient** — the vector of partial derivatives $\partial f/\partial x_j$ of a scalar function; the Jacobian of $f$, written as a column. It points uphill.
**Hessian** — the matrix of second partial derivatives $\partial^2 f/\partial x_i\partial x_j$; the Jacobian of the gradient. Symmetric.
**Positive definite matrix** — a symmetric matrix with all eigenvalues positive; a Hessian that is positive definite at a point with zero gradient marks a minimum.
**Saddle point** — a point with zero gradient where the Hessian has eigenvalues of both signs: uphill in some directions, downhill in others.
**Quadratic model** — the second-order Taylor approximation $f(\mathbf x + \mathbf s) \approx f + \nabla f^\mathsf T\mathbf s + \tfrac12\mathbf s^\mathsf TH\mathbf s$.
**Steepest descent** — iteration $\mathbf x_{k+1} = \mathbf x_k - \alpha_k\nabla f(\mathbf x_k)$; converges linearly at a rate set by the condition number of the Hessian.
**Step size (learning rate)** — the multiplier $\alpha$ of the search direction in one iteration.
**Line search** — choosing the step size along a given search direction.
**Backtracking (Armijo) line search** — start with $\alpha = 1$ and halve until $f$ decreases by at least $c\,\alpha$ times the predicted decrease (sufficient decrease, $c \approx 10^{-4}$).
**Descent direction** — a direction $\mathbf s$ with $\nabla f^\mathsf T\mathbf s < 0$, along which $f$ decreases for small steps.
**Zig-zagging** — the back-and-forth path of steepest descent across a long narrow valley (badly conditioned Hessian).
**Newton's method for minimisation** — solve $H\mathbf s = -\nabla f$ and step to the minimum of the quadratic model; quadratic convergence near a minimum.
**Modified Newton** — Newton with the Hessian shifted to $H + \mu I$ so that it is positive definite and every step goes downhill.
**Trust region** — a safeguard that limits each step to a region where the quadratic model is trusted, instead of a line search.
**Quasi-Newton method** — a method that builds an approximation of the Hessian (or its inverse) from changes of the gradient, without second derivatives.
**Secant condition** — the requirement $B\mathbf y = \mathbf s$ that the curvature model reproduces the observed change of gradient $\mathbf y$ along the last step $\mathbf s$.
**BFGS** — the standard quasi-Newton update (Broyden, Fletcher, Goldfarb, Shanno), a rank-2 update of the approximate inverse Hessian that keeps it positive definite.
**L-BFGS** — limited-memory BFGS: stores only the last few step and gradient-change pairs, for problems with very many parameters.
**Nelder–Mead method** — derivative-free minimisation that moves a simplex of $n + 1$ points by reflection, expansion, contraction and shrinking.
**Simplex** — $n + 1$ points in $n$ dimensions (a triangle in 2-D), as used by Nelder–Mead.
**Conjugate gradient (nonlinear)** — a first-order method that combines the new gradient with the previous direction, using only $O(n)$ memory.
**Parameter scaling** — rescaling parameters (units, logs, `parscale`) so that they are of similar size; it lowers the condition number of the Hessian.
**Multistart** — running a local optimiser from many starting values and comparing the results.
**Four-parameter logistic (4PL)** — the dose–response model $b + (t-b)/(1 + (\mathrm{EC50}/c)^h)$ with bottom, top, EC50 and Hill slope.
**Stochastic gradient descent (SGD)** — steepest descent using the gradient of a small random subset (minibatch) of the data terms at each step.

## Install

- No new packages. R: numDeriv (already listed by earlier units). Python: numpy, scipy. Rust: argmin, argmin-math, nalgebra (already in Cargo.toml).
