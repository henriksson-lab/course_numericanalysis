**Runge–Kutta method** — A one-step ODE method that evaluates the slope at several points inside the step and combines them with fixed weights.
**Stage** — One evaluation of the right-hand side f inside a Runge–Kutta step; the cost of a step is its number of stages.
**Midpoint method** — Second-order Runge–Kutta method that uses the slope at a half Euler step for the whole step.
**Heun's method** — Second-order Runge–Kutta method that averages the slopes at the start and at the Euler end point (explicit trapezoid rule).
**RK4 (classical Runge–Kutta)** — Fourth-order method with four slopes (start, midpoint twice, end) and weights 1, 2, 2, 1 over 6.
**Function evaluation** — One call of the ODE right-hand side f; the natural measure of the cost of an ODE solver.
**Embedded pair** — Two Runge–Kutta formulas of different order sharing the same stages; their difference estimates the local error.
**Local extrapolation** — Continuing with the higher-order answer of an embedded pair although the error estimate refers to the lower-order one.
**Adaptive step size** — Choosing each step size from an error estimate so that the estimated error per step meets a tolerance.
**Step rejection** — Discarding a step whose estimated error exceeds the tolerance and retrying it with a smaller step.
**Step-size controller** — The rule h_new = h · min(facmax, max(facmin, 0.9 · err^(-1/(q+1)))) that sets the next step from the scaled error.
**Dormand–Prince 5(4)** — The standard 7-stage embedded pair of orders 5 and 4 behind ode45, RK45 and Dopri5.
**FSAL (first same as last)** — Property of a pair whose last stage equals the first stage of the next step, saving one evaluation per step.
**Relative tolerance (rtol)** — The allowed error per step as a fraction of the size of each variable.
**Absolute tolerance (atol)** — The allowed error per step that applies when a variable is small; it has the units of the variable.
**Dense output** — A polynomial through each step, built from its stages, that gives the solution at any time inside the step.
**Event** — A time at which a chosen function g(t, y) of the solution crosses zero, located by the solver (e.g. an epidemic peak or a threshold).
**Conserved quantity (invariant)** — A function of the state that stays constant along exact solutions; its drift measures the error of a numerical solution.
**Lotka–Volterra model** — Predator–prey ODE model whose solutions are closed orbits around an equilibrium.
**Repressilator** — A synthetic ring of three repressor genes (Elowitz & Leibler 2000) that oscillates.

## Install

Nothing beyond intro/install.qmd (deSolve, numpy, scipy, matplotlib; Rust crate ode_solvers already in Cargo.toml).
