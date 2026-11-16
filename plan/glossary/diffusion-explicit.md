**Partial differential equation (PDE)** — Equation linking partial derivatives of a function of several variables, e.g. c(x, t) in space and time.
**Partial derivative** — Slope with respect to one variable while the others are held fixed; ∂c/∂x is the gradient at a fixed time, ∂c/∂t the rate of change at a fixed place.
**Fick's law** — Diffusive flux is proportional to the concentration gradient and runs downhill: J = −D ∂c/∂x.
**Diffusion coefficient** — The constant D in Fick's law, in m²/s; about 1e-9 m²/s (1000 µm²/s) for a small molecule and 1e-10 m²/s for a small protein in water.
**Diffusion equation (heat equation)** — ∂c/∂t = D ∂²c/∂x², from Fick's law and a mass balance; heat conduction obeys the same equation with thermal diffusivity α.
**Thermal diffusivity** — The heat-equation analogue of D, α = k/(ρ c_p); about 1.4e-7 m²/s for water.
**Diffusion time scale** — L²/D, the time for diffusion over a distance L; doubling L quadruples it.
**Initial condition** — The profile c(x, 0) at the start of a PDE simulation.
**Boundary condition** — What holds at the ends of the domain for all times; needed at each end of a 1D diffusion problem.
**Dirichlet boundary condition** — A fixed value at the boundary, e.g. a large reservoir holding the concentration constant.
**No-flux (Neumann) boundary condition** — Zero flux through the boundary, ∂c/∂x = 0, as at a wall.
**Ghost point** — A fictitious grid point outside the domain, set so that a boundary condition holds (c_{N+1} = c_{N−1} for no flux), letting the central stencil be used at the boundary.
**Grid (mesh)** — The points x_i = iΔx and times t_n = nΔt at which a PDE solution is approximated by values c_i^n.
**Method of lines** — Discretising a PDE in space only, which turns it into a large system of ODEs in time.
**FTCS scheme** — Forward-time, centred-space: Euler in time with the (1, −2, 1) stencil in space, c_i^{n+1} = c_i^n + r(c_{i−1}^n − 2c_i^n + c_{i+1}^n).
**Explicit scheme** — A time-stepping scheme in which each new value is computed directly from old values, without solving a system.
**Mesh ratio** — r = D Δt / Δx², the single number that governs the stability of FTCS.
**Stability limit of FTCS** — r ≤ 1/2, i.e. Δt ≤ Δx²/(2D); larger steps make a grid-scale zigzag grow without bound.
**Zigzag (sawtooth) mode** — The shortest wave on a grid, c_i = (−1)^i; the fastest-decaying mode of the diffusion grid and the one that limits explicit steps.
**Maximum principle** — Diffusion never creates values above the largest or below the smallest initial/boundary value; violating it signals instability or a bug.
**Second-difference matrix** — Tridiagonal matrix with −2 on the diagonal and 1 beside it; its eigenvalues −4 sin²(kπ/2N) give the decay rates of the grid modes.
**Nondimensionalisation** — Rescaling variables (x/L, tD/L²) so that the parameters drop out; all diffusion chambers become one problem.
**Fourier-series solution** — Exact solution of the diffusion equation as a sum of sine modes, each decaying at rate n²π²D/L².
**Grid refinement** — Repeating a simulation on finer grids to check that the answer converges at the expected order (here error ÷4 per halving of Δx).
**Chemotaxis chamber (diffusion-based gradient generator)** — A flow-free channel between a source and a sink reservoir in which diffusion sets up a linear concentration gradient for migrating cells.

## Install

None beyond `intro/install.qmd` (R: base/stats/graphics only; Python: numpy, matplotlib; Rust: standard library, plus `nalgebra` (already in Cargo.toml) for the eigenvalues in C4).
