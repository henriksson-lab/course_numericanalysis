# Audience

Biologists who want to understand the numerical machinery under the tools they
already use. The course is **self-paced with no time limit**: there is no
calendar to fit, so units are sized by topic, and the course can be as long as
the material needs.

## What we can expect them to know

### Biology / domain
- Molecular, cell and population biology; have seen real data (expression
  tables, time courses, dose–response, growth curves, images).

### Programming
- **Can program in R or Python** (not necessarily both): vectors/arrays, data
  frames, functions, loops, plotting, installing packages.
- Have *called* numerical routines (`lm`, `optim`, `nls`, `scipy.optimize`,
  `solve`) without knowing what happens inside.
- Caveat: many know only one of the two languages. Teaching is in **R**, so R
  idioms (vectorisation, `%*%`, 1-indexing, `function`) are explained the
  first time they matter. Slides are R only; exercise
  pages show every code example, hint and solution in R, Python and Rust.

### Mathematics
- **A first calculus course:** derivatives and integrals in one variable,
  the chain rule, Taylor polynomials, limits. Probably weak on: Taylor's theorem
  with remainder, multivariate calculus (partial derivatives, gradient, Jacobian,
  Hessian), and series convergence.
- **A first linear algebra course:** matrix product, inverse, $A\mathbf{x} =
  \mathbf{b}$ by Gaussian elimination, determinant, rank, eigenvalues of small
  matrices, maybe orthogonality. Probably weak on: norms other than Euclidean,
  factorisations as *objects* (LU, QR, SVD), and anything about cost or accuracy.
- **Not assumed:** real analysis, proofs, $\varepsilon$–$\delta$ arguments,
  complex numbers beyond the basics, differential equations beyond
  $y' = ky$, or any prior numerical analysis.

## What they have not seen

This is the gap the course closes:

- How numbers are **stored** (floating point) and what that costs.
- **Error** as a quantity to estimate: truncation vs rounding, conditioning vs
  stability, order of convergence.
- **Cost**: operation counts, $O(n^3)$ vs $O(n)$, why sparsity matters.
- What the black boxes do: iterative root finding, factorisations, adaptive
  step-size control, quasi-Newton optimisation — and how each one **fails**.

## Pedagogical implications
- Lead every method with a biological problem that needs it (see
  [designprinciples.md](designprinciples.md)).
- Recap calculus and linear algebra **just in time**, with numbers, and link a
  refresher. Multivariate calculus (gradient, Jacobian) is introduced in the
  course, where first needed (Newton for systems / optimisation).
- Prefer **experiments over proofs**: show the convergence order on a log–log
  plot rather than proving it. Derive only what fits on a slide and pays off.
- Make failure visible: every method gets a slide on how it breaks.
