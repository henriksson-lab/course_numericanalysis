# Brief for writing a unit (teacher-only, not rendered)

Use this when writing a new unit. Read first:

1. `designprinciples.md`, `audience.md`, and your unit's entry in `contentplan.md`
   (for Part 10 also the "Literature sources for Part 10" list there).
2. `plan/units.md` — canonical unit numbers, titles and parts. Use these exact
   titles whenever you mention another unit.
3. **Our own finished units** are the primary format reference — read
   `floating-point/lec.qmd`, `floating-point/exercise.qmd` (incl. its Rust
   tabs and `rust/examples/floating_point_*.rs`) and `floating-point/quiz.qmd`.
   Match their structure, tone, density and conventions.
4. Secondary reference implementation of the format (sister course, read-only — never edit it):
   `/data/henriksson/github/course/quantbio-3MB00A/eigen/lec.qmd`,
   `.../eigen/exercise.qmd`, `.../eigen/quiz.qmd`,
   `.../solving-systems-r/lec.qmd`, `.../compensation/make-comp-figs.R`.
   Copy the YAML headers (revealjs settings, `theme: [simple, ../theme.scss]`,
   1280×720, katex), the "About this lecture" opener, `::: notes` on every
   slide, `{.smaller}` where needed, and the exercise/quiz page structure.

## Files to produce, all inside `<slug>/`

- `lec.qmd` — reveal.js deck, **≤ 30 slides including the title slide**
  (aim 20–25). Ends with a Summary slide. Subtitle: `"Part N — <part name>"`,
  author: `"Numerical analysis for biologists"`.
  **Slides use R only** (plain ```` ```r ```` fences, outputs shown as
  comments `# [1] ...` or in a separate output fence). Each code block ≤ ~12 lines.
- `exercise.qmd` — Goal → Background → By hand (B1, B2, …) → In code (C1, C2, …)
  → Where this leads. Solutions and hints in collapsible callouts as in the
  reference. **Every code block on this page — background examples, prompts,
  hints, solutions — is an R + Python + Rust tabset**:

  ````md
  ::: {.panel-tabset group="language"}

  ### R

  ```r
  ...
  ```

  ### Python

  ```python
  ...
  ```

  ### Rust

  ```rust
  ...
  ```

  :::
  ````

  Python uses numpy/scipy/matplotlib only. Same printed output in both where
  possible. Rust: each snippet is a runnable file
  `rust/examples/<unit_snake>_<task>.rs` (see `rust/README.md` for naming and
  the allowed crates); the tab shows its code and real output. Snippets that
  are fragments in R/Python may be complete `fn main()` programs in Rust.
  Rust output goes in a trailing `// Output:` comment block in the tab (the file
  stays runnable); the first line of each Rust tab is
  `// rust/examples/<name>.rs  —  cargo run --example <name>` (add `--release`
  for timing). Write and run the Rust yourself in this pass.
  **No interactive demo / Explore section in this pass** (a later pass adds it).
- `quiz.qmd` — 8–12 short questions, each with a collapsible Hint and Answer.
  Code in quiz answers, if any, is also an R + Python + Rust tabset.
- `make-<slug>-figs.R` — generates every figure as SVG (`svg()` from base R
  grDevices; cairo is available, svglite and ggplot2 are not guaranteed —
  use base graphics). Figures named `figNN-short-name.svg` in the unit folder.
  Header comment like the reference: what it makes, how to run
  (`cd <slug> && Rscript make-<slug>-figs.R`). Labels large enough for slides
  (cex ≥ 1.3), colour-blind-safe colours.

## Code correctness

Code is **static** (not executed at render time), so you must **run every R
, every Python and every Rust snippet yourself** (`Rscript --vanilla`,
`python3 -I`, `cd rust && cargo run -q --example <name>`) and
paste the real output. Put scratch scripts in a scratch directory, not in the
unit folder. Floating-point outputs must be exactly what R/numpy print.
Available: R 4.6 with base/stats plus Matrix, irlba, RSpectra, minpack.lm,
numDeriv, deSolve, rootSolve, lpSolve, nnls, pracma, splines (no ggplot2 for
figures: base graphics); Python 3 with numpy 2.4, scipy 1.15, matplotlib 3.10;
Rust crates listed in `rust/README.md` (do not edit Cargo.toml — report a
missing crate instead, and hand-write the method meanwhile). Other workers
build in the same Cargo project concurrently; a build lock wait is normal.

## Content rules

- Pair every method with biology (see the unit entry). The biology must *need*
  the numerics.
- Show the failure: every method or concept gets a case where it goes wrong,
  and the check a student runs to detect it.
- By hand for small cases (exercise "By hand" section is pen and paper only).
- Notation: introduce slowly, numbers next to every formula.
- Language: plain, concise English; British spelling (optimise, behaviour);
  no hype words.
- Link the first mention of a key R function to its documentation
  (e.g. https://stat.ethz.ch/R-manual/R-devel/library/base/html/zMachine.html)
  and the Python equivalent in the exercise.
- Only write inside your own unit folder, `rust/examples/<unit_snake>_*`, and
  the two sidecar files below. Do not edit `_quarto.yml`, `dictionary.qmd`,
  `intro/install.qmd`, `contentplan.md` or other units.
- **Cross-references (two-pass approach).** In this pass, refer to other units
  by their exact title in italics (*Newton and secant methods*), optionally
  with the part (Part 2) — no hyperlinks. Assume earlier units (lower number)
  exist and their content is as in `contentplan.md`; you may read them if they
  exist. Do not re-teach their content beyond a one-slide recap. Record every
  cross-reference in `plan/xref/<slug>.md` as lines
  `- <file>:<line or section> -> <target slug> (<what is referred to>)`.
  A later pass turns these into links and checks consistency.
- **Shared running examples** (keep numbers consistent across units): logistic
  OD growth curve K = 1.2, N0 = 0.05, r = 0.8 per hour (see
  `convergence-cost`); if you introduce a new running example that later units
  should reuse, describe it in your xref file under a heading
  `## Running examples defined here`.
- **Further-reading units** (Part 11): add `<slug>/_metadata.yml` containing
  `include-before-body: ../_further-reading-banner.html` (as in the sister
  course), and they may be shorter (15–20 slides).

## Check before reporting

- Never kill processes you did not start (no `pkill`/`killall`): other workers run Chrome and builds too.

- `quarto render <slug>/lec.qmd`, `quarto render <slug>/exercise.qmd`,
  `quarto render <slug>/quiz.qmd` from the project root all succeed
  (if a render collides with another one running, wait and retry).
- Count slides: `grep -c '^## ' <slug>/lec.qmd` (+1 title) ≤ 30.
- Check no slide is overfull: rough rule ≤ ~10 lines of bullets or one figure
  plus ≤ 4 bullets, or code ≤ 12 lines plus ≤ 2 bullets.

## Report back

Write `plan/glossary/<slug>.md` with glossary terms
(`**Term** — one-line definition.` per line) and, at the end, a section
`## Install` listing any R/Python packages beyond those already on
`intro/install.qmd`. Then reply briefly: files written; slide count; anything
you were unsure about or deferred (≤ 10 lines).
