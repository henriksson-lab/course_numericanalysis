# Design principles

The principles below govern how we build course material in this repository. They
exist to keep us honest as the material scales — every new exercise, lecture, and
chapter should be checked against them.

See [audience.md](audience.md) for who we are designing for and
[contentplan.md](contentplan.md) for the planned units.

## 1. Method ↔ biology pairing

Every mathematical or statistical method is introduced together with a biological
problem that motivates it. We do not teach a technique in the abstract first and
apply it to biology later.

- The unit of design is the **pairing**. There is no calendar: shift, merge or
  split units as the pairings become clearer.
- The biology must *need* the numerics (an underflowing likelihood, a stiff
  kinetic model, a fit that converges to the wrong minimum), not just decorate it.

## 2. Show the math *and* the code

Every formula is accompanied by code that computes it, and every code result is
tied back to the math it implements. A student should never see a derivation
without being able to run it, nor a function call without knowing what it computes.

- **R is the main language. Slides are R only** — we do not have the
  manpower to maintain Python versions of the decks.
- **Exercise pages also have Python and Rust.** Code on exercise and quiz
  pages (background examples, prompts, hints, solutions, quiz answers) is a
  Quarto `panel-tabset` with an **R** tab first (shown by default), then
  **Python** (numpy/scipy/matplotlib), then **Rust**. Prose and figures follow
  the R version. Python and Rust users read the R on the slides and do the
  exercises in their own language.
- **Rust code lives in `rust/examples/`** as one runnable file per snippet
  (`<unit>_<task>.rs`, see `rust/README.md`); the Rust tab shows that file and
  its real output. Use the standard library when it is enough, otherwise only
  the crates listed in `rust/README.md`. Where Rust has no counterpart to a
  library call (e.g. R's `integrate`), name the crate used, or say plainly
  that the Rust version is hand-written.
- Keep the two tabs equivalent line by line where possible and make them print
  the same output. Where the library call differs (`uniroot` vs
  `scipy.optimize.brentq`, `deSolve::ode` vs `solve_ivp`), name both, and link
  both in hints.
- When a slide introduces a library call, the matching exercise names the
  Python equivalent, so Python users can map the slide to their language.
- Both versions are run before they are published, not just written.
- Only one-line inline mentions in prose (`solve()`) may be R only.

  Template:

  ````md
  ::: {.panel-tabset group="language"}

  ### R

  ```r
  x <- seq(0, 1, length.out = 11)
  ```

  ### Python

  ```python
  x = np.linspace(0, 1, 11)
  ```

  ### Rust

  ```rust
  let x: Vec<f64> = (0..=10).map(|i| i as f64 / 10.0).collect();
  ```

  :::
  ````

  `group="language"` makes all tabsets on a page switch together, so a reader
  who picks Python or Rust once sees that language everywhere on the page.
- **Write it yourself first, then use the library.** The first encounter with a
  method is a short hand-written implementation (a 10-line Newton, a 15-line
  RK4); then compare with the production routine (`uniroot`/`brentq`,
  `deSolve::ode`/`solve_ivp`, `optim`/`minimize`) and explain what it does extra.
- Build geometric and probabilistic *intuition* (pictures, small worked numbers)
  before notation.

## 3. Exercises over stretch goals

We prefer a good number of core exercises that every student is expected to
complete over a few exercises with long optional "going further" tails.

- The audience is broad. Stretch goals systematically benefit students who are
  already ahead and widen the spread rather than building shared competence.
- A short, neutral **"Where this leads"** note at the end of an exercise is fine —
  it tells students what we defer and to which unit, without inviting them to
  attempt it now.

## 4. Linked hints over named hints

Every key function, concept, or distribution named in a hint links to its
documentation — the R help page, the package vignette, or a canonical reference —
so a student never has to guess search terms to finish a core exercise.

## 5. By hand for the small case, in code for the large

Core exercises build the ability to *do the small case with pen and paper* — two
Newton steps, an LU of a 3×3, one Euler step, a trapezoid sum with three points,
a float rounded to binary — so the student knows what the code is doing. Then
the same idea is run in code at a size where only code works.

## 5b. Make error and failure visible

Every method gets: (1) an **error-vs-step** or **error-vs-iteration** plot on log
axes, from which the order is read off; (2) at least one example where it
**fails** (diverges, stalls, is unstable, returns a confident wrong answer); and
(3) the check a student should run to detect that failure in their own work.

## 6. Honest about what we defer and what we assume

When an exercise uses a simplification we will later replace (a tiny dataset, a
hand-rolled solver instead of the optimised one, a fixed step size before adaptive
steps), say so plainly and name the unit in which the proper treatment arrives.
Students should always know what is a stepping stone and what is the real
thing.

## 7. Keep slide decks small — split, don't cut

Each lecture deck may be recorded as a **video**, so a deck should be a
self-contained snippet on one topic.

No single lecture deck should grow large. When new material would make a deck
too big, first ask whether it is really **its own topic** and should become a
**separate deck** (its own unit folder with `lec.qmd` and `exercise.qmd`, placed
next to the original in the sidebar). Prefer splitting over cutting good
material. Cut or merge only slides that genuinely repeat a point.

**Hard ceiling: 30 slides** per deck, title slide included. Aim for 20–25; the ceiling is not a target.

## 7b. Figures stay small

Line art is SVG; figures with many points or dense images (big scatters,
bifurcation diagrams, ensembles of trajectories, heat maps) are PNG at
150 dpi. No figure file over about 1 MB (GitHub rejects files over 100 MB,
and big SVGs make pages slow).

## 8. Slide and exercise format

Follow the conventions of the sister course (quantbio-3MB00A): every deck opens
with an "About this lecture" slide, every slide has `::: notes` speaker notes,
the exercise page has Goal → Background → By hand → In code → Explore (demo),
and the quiz has a collapsible hint and answer per question.
