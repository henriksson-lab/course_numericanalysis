**Partial derivative** — the derivative of a function of several variables with respect to one of them, all others held fixed; written $\partial F/\partial x_j$.
**Jacobian** — the matrix of all first partial derivatives of a vector function $\mathbf F$: entry $(i, j)$ is $\partial F_i/\partial x_j$ (row = equation, column = variable).
**Linearisation** — the approximation $\mathbf F(\mathbf x + \boldsymbol\delta) \approx \mathbf F(\mathbf x) + J(\mathbf x)\,\boldsymbol\delta$, valid for small steps $\boldsymbol\delta$.
**Newton's method for systems** — iteration that solves $J(\mathbf x_k)\boldsymbol\delta = -\mathbf F(\mathbf x_k)$ and sets $\mathbf x_{k+1} = \mathbf x_k + \boldsymbol\delta$; one linear solve per step.
**Residual** — the value $\mathbf F(\mathbf x)$ at an approximate root; its length $\|\mathbf F(\mathbf x)\|$ measures how far the equations are from being satisfied.
**Euclidean norm** — the length of a vector, $\|\mathbf a\| = \sqrt{a_1^2 + \dots + a_n^2}$.
**Singular matrix** — a square matrix with determinant 0, which has no inverse; a singular Jacobian is the system analogue of $f'(x) = 0$.
**Finite-difference Jacobian** — an approximate Jacobian built column by column as $(\mathbf F(\mathbf x + h_j\mathbf e_j) - \mathbf F(\mathbf x))/h_j$, with $h_j \approx \sqrt\varepsilon\,\max(|x_j|, 1)$.
**Steady state** — a state where all rates of change are zero, i.e. a root of the right-hand side of an ODE system.
**Bistability** — the existence of two stable steady states for the same parameters, so the system's state depends on its history (cellular memory).
**Toggle switch** — a two-gene circuit in which each gene's protein represses the other gene, giving bistability for strong enough promoters.
**Hill function** — a sigmoidal dose–response $\alpha/(1 + x^n)$ (repression) or $\alpha x^n/(1 + x^n)$ (activation); $n$ is the Hill coefficient (cooperativity).
**Basin of attraction (of a root finder)** — the set of starting points from which an iteration converges to a given root.
**Damped Newton** — Newton's method with step length $t\boldsymbol\delta$, $0 < t \le 1$, shortened until the step is acceptable (e.g. the residual does not grow).
**Trust region** — a safeguard that limits each step to a region where the linear model is trusted; used by MINPACK's hybrid method (`scipy.optimize.root`).

## Install

- R: `rootSolve` (multiroot), `numDeriv` (jacobian, used once as a comparison in the exercise).
- Python: none beyond numpy/scipy.
- Rust: none beyond `nalgebra`.
