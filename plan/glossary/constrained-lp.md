**Constraint** — A restriction on the parameters of an optimisation problem, either an equality ($\sum_i p_i = 1$) or an inequality ($c \ge 0$).
**Feasible set** — The set of parameter values that satisfy all constraints.
**Active constraint** — An inequality constraint that holds with equality at a solution; inactive constraints do not change the answer.
**Bound constraint** — A constraint of the form $\ell \le \theta \le u$ on a single parameter, handled by L-BFGS-B.
**L-BFGS-B** — Limited-memory BFGS with bound constraints: steps that would leave the box are projected back onto it.
**Reparameterisation** — Optimising a free parameter $\theta$ mapped into the allowed set, e.g. $c = e^\theta$ for positivity or softmax for proportions.
**Softmax** — The map $p_i = e^{\theta_i}/\sum_j e^{\theta_j}$, which turns any real vector into positive proportions summing to 1.
**Lagrange multiplier** — The number $\lambda$ in $\nabla f = \lambda \nabla g$ at a constrained optimum; the rate at which the optimal value changes when the constraint level changes.
**Shadow price** — The change in the optimal objective per unit relaxation of a constraint; zero for an inactive constraint.
**KKT conditions** — Optimality conditions for inequality-constrained problems: gradient balanced by multipliers, multipliers $\ge 0$, and zero for inactive constraints.
**Linear program (LP)** — An optimisation problem with a linear objective and linear equality and inequality constraints.
**Feasible polytope** — The convex polygon (in higher dimensions, polytope) of points satisfying the constraints of an LP.
**Vertex (corner)** — A point of the feasible polytope where as many independent constraints are active as there are variables; an LP optimum, if it exists, is attained at one.
**Simplex method** — Dantzig's LP algorithm that moves from vertex to vertex along edges that improve the objective until no edge does.
**Ratio test** — The simplex step that finds how far one can move along an edge before the first constraint becomes active.
**Interior-point method** — An LP (or nonlinear) algorithm that follows a path through the inside of the feasible set using Newton steps.
**Degenerate vertex** — A vertex where more constraints are active than needed to define it; shadow prices there are not unique.
**Unbounded LP** — An LP whose objective can be improved without limit, usually because a constraint is missing.
**Infeasible LP** — An LP whose constraints cannot all be satisfied at once.
**Stoichiometric matrix** — The matrix $S$ with one row per metabolite and one column per reaction, entry = molecules made (+) or used (−).
**Flux** — The rate of a reaction in a metabolic network (e.g. mmol per gram dry weight per hour).
**Steady state (FBA)** — The assumption $S\mathbf v = \mathbf 0$: every internal metabolite is produced as fast as it is consumed.
**Exchange reaction** — A pseudo-reaction for uptake from or secretion to the medium; its bounds encode nutrient availability.
**Biomass reaction** — A pseudo-reaction draining precursors and energy in the proportions needed to make new cells; its flux is the growth rate.
**Flux balance analysis (FBA)** — Maximising an objective (usually the biomass flux) subject to $S\mathbf v = \mathbf 0$ and flux bounds, as a linear program.
**Flux variability analysis (FVA)** — Fixing the objective at its optimum and minimising and maximising each flux, giving the range of each flux over all optimal solutions.
**Synthetic lethal pair** — Two reactions (genes) each dispensable alone but lethal when both are removed.
**Overflow metabolism** — Secretion of fermentation products (e.g. acetate) even when respiration is possible, here because oxygen limits respiration.
**COBRA / cobrapy** — Standard toolboxes (MATLAB / Python) for constraint-based modelling of genome-scale metabolic networks.

## Install

R: `lpSolve` (LP solver) and `numDeriv` (numerical gradient) — neither is on `intro/install.qmd` yet. Python: no new packages (scipy's `linprog` uses HiGHS). Rust: `minilp` and `argmin` (already in Cargo.toml).
