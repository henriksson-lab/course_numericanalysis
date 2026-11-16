**Ordinary differential equation (ODE)** — Equation for the rate of change of a quantity, y' = f(t, y), in terms of the quantity itself (and time).
**Initial value problem (IVP)** — An ODE together with a start value y(t0) = y0; its solution is the curve that starts there and follows the slope f.
**Rate in − rate out** — Modelling principle: write the derivative of each quantity as the sum of flows into it minus flows out of it.
**Steady state (equilibrium)** — A value y* where f(y*) = 0, so the solution stays constant if it starts there.
**Direction field (slope field)** — Picture of an ODE: a short line with slope f(t, y) at each point of the (t, y) plane; solutions are tangent to it.
**Autonomous ODE** — ODE whose right-hand side does not depend on t explicitly, y' = f(y).
**Logistic growth** — N' = rN(1 - N/K): growth at per-capita rate r that slows to zero at the carrying capacity K.
**Carrying capacity** — The stable equilibrium K of the logistic model, the density at which growth stops.
**Euler's method** — Explicit one-step ODE method y_{n+1} = y_n + h f(t_n, y_n): follow the tangent for a step h.
**Step size** — The time increment h between successive approximations y_n ≈ y(t_n).
**Local error** — Error made in a single step started on the exact solution; O(h^2) for Euler.
**Global error** — Difference y_n - y(t_n) after all steps from the start; O(h) for Euler.
**Order of a method** — The power p with global error ≈ C h^p; the slope on a log–log plot of error against h. Euler has order 1.
**Step halving** — Running with h and h/2 and using the difference as an error estimate (for order 1, y_h - y_{h/2} ≈ error of y_{h/2}).
**Numerical instability (of an ODE method)** — Growth of errors from step to step because the step is too large, e.g. |1 - kh| > 1 for Euler on y' = -ky.
**Stability bound** — The largest stable step, h < 2/k for Euler on y' = -ky.
**Linearisation** — Replacing f near an equilibrium y* by f'(y*)(y - y*); gives the local decay rate k = -f'(y*) that sets the stability bound.
**Spurious oscillation / cycle** — An oscillation or cycle produced by the numerical method, absent from the true solution; it disappears when the step is reduced.
**System of ODEs** — Several coupled ODEs written as one vector equation y' = f(t, y); Euler applies unchanged with vectors.
**SIR model** — Epidemic model with susceptible, infected and recovered compartments: S' = -βSI/N, I' = βSI/N - γI, R' = γI.
**Basic reproduction number R0** — Expected number of new infections caused by one infected person in a fully susceptible population; β/γ in the SIR model.
**Conserved quantity (invariant)** — A combination of the variables that the exact solution keeps constant (S + I + R; an oscillator's amplitude); useful as a check.
**Higher-order ODE** — ODE involving second or higher derivatives; rewritten as a first-order system by naming x, x', x'', … as new variables.
**Harmonic oscillator** — x'' = -ω²x, solutions cos and sin with period 2π/ω; Euler makes its amplitude grow for every step size.
**Phase plane** — Plot of one variable against another (e.g. x against x'/ω) instead of against time.
**Stiffness** — Situation where a fast rate forces a small step for stability long after the fast process is over (treated in a later unit).

## Install

None beyond `intro/install.qmd` (R: base/stats/graphics only; Python: numpy, matplotlib; Rust: standard library only, no crates).
