**Reaction–diffusion equation** — a PDE in which a concentration both diffuses and is produced or consumed locally, $\partial c/\partial t = D\nabla^2 c + f(c)$.
**Synthesis–diffusion–degradation (SDD) model** — a morphogen made at a source, spreading by diffusion and degraded at rate $k$; gives an exponential steady gradient.
**Decay length** — $\lambda = \sqrt{D/k}$, the distance over which an SDD gradient falls by a factor $e$.
**Morphogen** — a signalling molecule whose concentration tells cells their position (e.g. Bicoid, Dorsal, Fgf8).
**Local accumulation time** — $\tau(x) = \int_0^\infty (1 - c(x,t)/c_s(x))\,dt$, a mean time for the gradient at $x$ to reach steady state.
**Ghost point** — a fictitious grid point outside the domain used to impose a flux or no-flux boundary condition to second order.
**Mass balance (gradient)** — at steady state total degradation $k\int c\,dx$ equals the influx $j$; a check that needs no exact solution.
**IMEX scheme** — implicit–explicit time stepping: stiff linear terms (diffusion) implicit, non-stiff or nonlinear terms (reaction) explicit.
**Operator splitting** — advancing each part of an equation (diffusion, reaction) separately within a time step.
**Lie splitting** — one full step of each part in turn; first-order accurate.
**Strang splitting** — half step of one part, full step of the other, half step of the first; second-order accurate.
**Splitting error** — the error of a splitting scheme even with exact sub-steps; zero only when the parts commute.
**Commutator** — $AB - BA$; zero when the order of applying two operators does not matter.
**Turing instability** — a uniform steady state that is stable without diffusion becomes unstable with diffusion at some wavelengths.
**Activator–inhibitor (activator–substrate) system** — two species where one promotes its own production and a faster-diffusing partner limits it.
**Schnakenberg model** — a simple two-species reaction ($a - u + u^2v$, $b - u^2v$) that shows Turing patterns.
**Wavenumber** — $q = 2\pi/\text{wavelength}$; a cosine perturbation $\cos(qx)$ turns $\nabla^2$ into $-q^2$.
**Dispersion relation** — the growth rate (largest real eigenvalue part) of a perturbation as a function of wavenumber.
**Numerical artefact (pattern)** — structure in a simulation produced by the grid or time step rather than the model; it changes when $\Delta x$ or $\Delta t$ is refined.
**Periodic boundary condition** — the domain wraps around (left edge joins right edge); natural for rings such as an embryo's circumference.

## Install

No packages beyond `intro/install.qmd` (uses R `Matrix::expm`, Python `scipy.linalg.expm`, Rust `nalgebra` and `rustfft`).
