**Finite volume method** — Discretisation that writes one conservation balance (flux in − flux out = source) per cell of a mesh; conserves mass exactly.
**Finite element method (FEM)** — Discretisation that writes the unknown as a sum of simple basis functions on a mesh and asks the equation to hold in a weighted-average (weak) sense.
**Mesh** — A division of the domain into elements (intervals, triangles, tetrahedra), possibly of different sizes, that can follow curved walls.
**Element** — One piece of a mesh (an interval in 1D, a triangle in 2D) on which the approximation is a simple polynomial.
**Node** — A mesh point where an unknown value lives; in linear elements, the corners of the elements.
**Hat function** — Piecewise linear basis function equal to 1 at its own node and 0 at all others.
**Weak form** — The equation multiplied by a test function and integrated (with integration by parts), so that only first derivatives appear.
**Natural boundary condition** — A flux condition (e.g. no flux) that the weak form satisfies automatically, with no extra equation.
**Element matrix** — The small matrix (2 × 2 in 1D, 3 × 3 for linear triangles) that one element contributes to the global system.
**Assembly** — Adding all element matrices into the rows and columns of their nodes to build the global sparse matrix.
**Stiffness matrix** — The assembled matrix of the diffusion term, $\int D\varphi_i'\varphi_j'$.
**Mass lumping** — Replacing the consistent reaction (mass) matrix by a diagonal one; in 1D on a uniform mesh it gives back finite differences.
**Graded mesh** — A mesh with small elements where the solution changes fast and large elements elsewhere.
**Adaptive mesh refinement** — Automatically splitting the elements where an error estimate is largest and solving again.
**Reynolds number** — $Re = \rho U \ell/\mu$, the ratio of inertial to viscous forces in a flow; small in microfluidics.
**Stokes flow** — Flow at $Re \ll 1$, where inertia is negligible: $\mu\nabla^2\mathbf u = \nabla p$, $\nabla\cdot\mathbf u = 0$; linear and laminar.
**Navier–Stokes equations** — The equations of mass and momentum conservation for a viscous fluid, including the nonlinear inertia term.
**Poiseuille flow** — Pressure-driven flow in a straight channel with a parabolic velocity profile, peak 1.5 × mean between plates.
**COMSOL Multiphysics** — Commercial graphical finite element package, widely used for microfluidic chip design.
**FEniCS / Firedrake** — Open-source Python finite element frameworks in which the weak form is written almost as on paper.

## Install

None beyond `intro/install.qmd` (base R; numpy/scipy/matplotlib with `scipy.linalg.solve_banded`; Rust std only). COMSOL, FEniCS, Firedrake, scikit-fem and deal.II are mentioned, not used.
