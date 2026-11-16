**Phase plane** — the plane of the two state variables of a two-variable ODE, on which the flow is drawn as arrows and trajectories.
**Trajectory (orbit)** — the curve traced in state space by one solution of an ODE as time runs.
**Vector field (direction field)** — the velocity (right-hand side of the ODE) drawn as an arrow at each point of state space.
**Nullcline** — the curve on which one variable does not change (e.g. du/dt = 0); the flow crosses it vertically or horizontally.
**Fixed point (equilibrium, steady state)** — a state where all derivatives are zero, F(x*) = 0.
**Stable fixed point** — small perturbations die out; for an ODE, all eigenvalues of J have negative real part.
**Linearisation** — replacing the ODE near a fixed point by dδ/dt = J δ, with J the Jacobian at the fixed point.
**Node** — fixed point with real eigenvalues of the same sign; trajectories approach (stable) or leave (unstable) without rotating.
**Saddle** — fixed point with real eigenvalues of opposite sign; attracting in one direction, repelling in another; unstable.
**Spiral (focus)** — fixed point with complex eigenvalues; trajectories rotate while approaching (Re λ < 0) or leaving (Re λ > 0).
**Centre** — fixed point with purely imaginary eigenvalues; linearisation predicts closed orbits but cannot decide stability.
**Trace–determinant plane** — classification of 2×2 fixed points by τ = tr J and Δ = det J: Δ < 0 saddle, τ < 0 < Δ stable, τ² < 4Δ spiral.
**Basin of attraction** — the set of starting states from which the ODE ends at a given stable state.
**Separatrix** — the curve (stable manifold of a saddle) separating two basins of attraction; computed by integrating backwards in time from the saddle.
**Stable manifold** — the set of states whose trajectories converge to a given saddle.
**Backward integration** — solving an ODE with decreasing time (equivalently, dx/dt = −F forwards); turns attracting directions into repelling ones.
**Limit cycle** — an isolated closed trajectory; nearby trajectories approach it (stable) and the system oscillates with a fixed period and amplitude.
**Bistability** — two stable steady states for the same parameters; the state reached depends on history.
**Bifurcation** — a parameter value at which the number or stability of steady states (or cycles) changes; an eigenvalue crosses Re λ = 0 there.
**Pitchfork bifurcation** — one stable state becomes unstable and two new stable states branch off (the symmetric toggle at α = 2).
**Hopf bifurcation** — a pair of complex eigenvalues crosses the imaginary axis and a limit cycle is born or dies.
**Bifurcation diagram** — steady states (and cycles) plotted against a parameter, stable and unstable branches drawn differently.
**Numerical continuation** — following a steady state as a parameter changes: predict from the previous solution, correct with Newton, monitor eigenvalues.
**Critical slowing down** — near a bifurcation the slowest eigenvalue approaches 0, so perturbations decay very slowly and simulations take very long to settle.
**Rightmost eigenvalue** — the eigenvalue with the largest real part; it decides stability (not the eigenvalue of largest modulus).

## Install

- R: `rootSolve` (`multiroot`, `jacobian.full`) is not listed on `intro/install.qmd` (also used by *Newton for systems of equations*). `RSpectra` is only mentioned on a slide, not run.
- Python and Rust: nothing new (numpy/scipy/matplotlib; `ode_solvers` and `nalgebra` are already in Cargo.toml).
