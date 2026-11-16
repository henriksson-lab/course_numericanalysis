**Steady state** — a solution that no longer changes in time; for diffusion, set dc/dt = 0.
**Boundary value problem (BVP)** — a differential equation with conditions on the boundary of the domain instead of an initial condition.
**Poisson equation** — D times the Laplacian of c equals a given source or sink term.
**Laplace equation** — the Poisson equation without sources: the Laplacian of c is zero.
**Laplacian** — the sum of the second derivatives in all directions, d2c/dx2 + d2c/dy2.
**Five-point stencil (five-point Laplacian)** — the finite-difference Laplacian on a 2D grid: (sum of the four neighbours - 4 times the centre) / h^2.
**Kronecker product** — A (x) B: the block matrix with blocks a_pq B; builds 2D grid matrices from 1D ones.
**Column-major / row-major** — storing a matrix column by column (R, Fortran) or row by row (numpy default, C); decides which index runs fastest when a grid is flattened.
**Manufactured solution** — a chosen smooth function whose derivatives are known exactly, used to test a discretisation.
**Maximum principle** — without sources, a diffusion solution cannot exceed its largest or fall below its smallest boundary value.
**Zero-order uptake** — consumption at a constant rate, independent of concentration.
**Michaelis–Menten uptake** — consumption R0 c / (Km + c): constant at high c, proportional to c at low c.
**Spurious root** — a solution of the discrete equations that has no physical meaning (e.g. negative concentration).
**Clipping (projected Newton step)** — replacing values outside the physical range (e.g. c < 0) by the boundary of that range after each step.
**Fill-reducing ordering** — a renumbering of the unknowns (AMD, nested dissection) that reduces fill-in in a sparse factorisation.
**Mother machine** — a microfluidic device with dead-end trenches, each holding a line of bacteria fed from a channel at the open end.
**Hypoxia** — oxygen below the level cells need for normal function (here: below 10 µM, about 1 % O2).

## Install

None. Matrix ships with R; scipy.sparse is part of scipy; sprs is in rust/Cargo.toml (it has no direct sparse solver, so the Rust solutions use a hand-written CG / Thomas).
