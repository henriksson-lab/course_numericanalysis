# Cross-references from `constrained-lp` (Constraints and linear programming)

- [ok: lu numbers checked; linked] constrained-lp/lec.qmd:46 ("Three kinds of constraint") -> lu (isotope mixing model, consumer (-24, 9) gives diet (0.55, -0.20, 0.65))
- [ok: svd -0.09 checked; linked] constrained-lp/lec.qmd:47 ("Three kinds of constraint") -> svd (least squares gives CD8 proportion -0.09)
- [fixed: notes add that these decay data differ from optimisation-1d's mRNA decay example] constrained-lp/lec.qmd:72 (notes, "A fit that wants a negative background") -> gradient-methods (optim BFGS)
- [ok] constrained-lp/lec.qmd:89 (notes, "Fix 1: bounds") -> gradient-methods (BFGS / L-BFGS)
- [fixed: identifiability does not discuss bounds; slide now says its SEs assume an interior minimum; linked] constrained-lp/lec.qmd:101 ("At a bound, the gradient is not zero") -> identifiability (Hessian standard errors invalid at a bound)
- [ok] constrained-lp/lec.qmd:104 (notes, same slide) -> svd (NNLS, KKT), identifiability (Hessian-based SEs)
- [ok: linked] constrained-lp/lec.qmd:111 ("Fix 2: reparameterise") -> nonlinear-least-squares (log-parameterisation of positive rates)
- [ok] constrained-lp/lec.qmd:166 ("The multiplier is a shadow price") -> svd (NNLS as KKT example)
- [ok: linked] constrained-lp/lec.qmd:232 ("How many corners?") -> newton-systems (interior-point methods take Newton steps)
- [ok] constrained-lp/lec.qmd:224 (notes, "The simplex idea") -> lu (corner = solve of active constraints; LU updates)
- [fixed: lu has no null space; slide says elimination (as in lu) leaves 3 free fluxes; linked; solution C5 likewise] constrained-lp/lec.qmd:302 ("The stoichiometric matrix") -> lu (null space of S by elimination)
- [ok: linked] constrained-lp/lec.qmd:343 ("The optimum is not unique: flux variability") -> identifiability (non-identifiability analogy)
- [ok: linked] constrained-lp/lec.qmd:411 ("Production tools") -> structured-sparse (sparse S)
- [ok: linked] constrained-lp/lec.qmd:426 (Summary) -> diffusion-explicit (next unit, Part 10)
- [ok: linked] constrained-lp/exercise.qmd:124 (B3) -> lu (isotope mixing model, consumer (-24, 9))
- [fixed: reworded, see lec.qmd:302] constrained-lp/exercise.qmd:951 (Solution C5) -> lu (null space of S)
- [ok: linked] constrained-lp/exercise.qmd:1382 (Where this leads) -> svd (NNLS and KKT)
- [ok: linked] constrained-lp/exercise.qmd:1383 (Where this leads) -> identifiability (parameter at bound, FVA range)
- [ok: linked] constrained-lp/exercise.qmd:1384 (Where this leads) -> diffusion-explicit (Part 10 start; sparse systems)
- [note: Part 8 named, not linked] constrained-lp/exercise.qmd:1385 (Where this leads) -> euler / runge-kutta (Part 8 ODEs; dynamic FBA couples FBA to an ODE)
- [ok] constrained-lp/quiz.qmd:189 (Q10) -> identifiability (non-identifiable parameter)

## Running examples defined here

- **mRNA decay with background** (synthetic): t = 0:8 h, y = (1.020, 0.577, 0.378, 0.243, 0.125, 0.062, 0.020, 0.010, -0.012) (true A = 1, k = 0.5, c = 0 plus fixed noise). Model A e^{-kt} + c. Free BFGS fit (1.0366, 0.4788, -0.0294), SSE 0.0023; with c >= 0 (L-BFGS-B) (1.0152, 0.5214, 0), SSE 0.0038, dSSE/dc = 0.105 at the bound; half-life 1.448 h vs 1.329 h.
- **ATP LP (overflow metabolism caricature)**: max 6 x1 + 3 x2 (respiration, fermentation), x1 + x2 <= 6 (carbon), 4 x1 + x2 <= 12 (proteome); optimum (2, 4), ATP 24, shadow prices 2 and 1. With oxygen x1 <= 1.5: (1.5, 4.5), 22.5, prices (3, 0, 3).
- [ok] **Toy FBA network** (9 reactions, 6 metabolites G, O2, X, P, ATP, Ac): EX_glc (<= 10), EX_o2 (<= 15), GLY: G -> 2P + 2ATP, BYP1: G -> X, BYP2: X -> 2P + 2ATP, RES: P + 2 O2 -> 11 ATP, FER: P -> Ac + ATP, EX_ac: Ac ->, BIO: P + 5 ATP ->. Rank S = 6. Max growth 15 (RES 5, O2 uptake 10), glucose shadow price 1.5, O2 0; FVA: GLY, BYP1, BYP2 in [0, 10]; KO RES -> 6.667, GLY+BYP1 synthetic lethal; O2 limit 4 -> growth 10, acetate 8, prices 0.667 / 0.833; acetate stops at O2 = 10. Could be reused in a later unit on sparse LP or dynamic FBA.
