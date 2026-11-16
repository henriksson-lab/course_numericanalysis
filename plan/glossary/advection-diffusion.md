**Advection** — transport of a dissolved substance by the flow of the fluid that carries it.
**Advection–diffusion equation** — $\partial c/\partial t + u\,\partial c/\partial x = D\,\partial^2 c/\partial x^2$: mass balance with a flux $uc - D\,\partial c/\partial x$.
**Péclet number** — $\text{Pe} = uL/D$, the ratio of the diffusion time $L^2/D$ to the flow time $L/u$; Pe ≫ 1 means advection dominates.
**Cell Péclet number** — $u\Delta x/D$, the Péclet number of one grid cell; central differences wiggle when it exceeds 2.
**Mixing length** — distance $\approx \bar u w^2/D$ that side-by-side streams travel before diffusion mixes them across a channel of width $w$.
**Boundary layer** — thin region (here of thickness $D/u$) where a solution changes steeply; the grid must resolve it.
**Central difference (advection)** — $(c_{i+1} - c_{i-1})/(2\Delta x)$ for $\partial c/\partial x$; second order but can give oscillating solutions.
**Upwind difference** — one-sided difference taken on the side the flow comes from, e.g. $(c_i - c_{i-1})/\Delta x$ for $u > 0$; first order, never wiggles.
**Wiggles** — spurious point-to-point oscillations in a discrete solution, here from central differences with cell Péclet number > 2.
**Numerical diffusion** — artificial extra diffusion introduced by a scheme; upwind adds $u\Delta x/2$ (steady) or $u\Delta x(1 - C)/2$ (explicit in time).
**Courant number** — $C = u\Delta t/\Delta x$, the number of grid cells the flow moves in one time step.
**CFL condition** — Courant–Friedrichs–Lewy stability limit $C \le 1$ for explicit advection; with diffusion, $C + 2r \le 1$ for explicit upwind.
**Poiseuille profile** — the parabolic velocity profile $u(y) = 6\bar u (y/w)(1 - y/w)$ of slow pressure-driven flow between two walls.
**Plug flow** — idealised flow with the same speed across the whole channel.
**Marching** — solving a steady problem with Pe ≫ 1 step by step in the downstream coordinate, treating it like time.
**T-sensor** — microfluidic device in which two streams flow side by side; the diffusion interface between them is used to measure concentrations or diffusion coefficients.
**H-filter** — membrane-free microfluidic filter: small molecules diffuse from a sample stream into a buffer stream while cells and large molecules stay behind.
**Christmas-tree gradient generator** — microfluidic network that repeatedly splits, mixes and recombines inlet streams to produce a concentration gradient of a chosen shape.

## Install

None beyond `intro/install.qmd` (R: base/stats; Python: numpy, scipy; Rust: `statrs` and `argmin`, already in `Cargo.toml`).
