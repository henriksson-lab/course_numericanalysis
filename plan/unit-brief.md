# Brief for writing a unit (teacher-only, not rendered)

Use this when writing a new unit. Read first:

1. `designprinciples.md`, `audience.md`, and your unit's entry in `contentplan.md`.
2. Reference implementation of the format (sister course, read-only — never edit it):
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
  Interactive demo: see `plan/demo-brief.md`.
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
Available: R 4.6 (base + deSolve may not be installed: use base R only),
Python 3 with numpy 2.4, scipy 1.15, matplotlib 3.10.

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
- Only write inside your own unit folder. Do not edit `_quarto.yml`,
  `dictionary.qmd`, `intro/install.qmd` or other units; instead, report
  what should be added there.

## Check before reporting

- `quarto render <slug>/lec.qmd`, `quarto render <slug>/exercise.qmd`,
  `quarto render <slug>/quiz.qmd` from the project root all succeed
  (if a render collides with another one running, wait and retry).
- Count slides: `grep -c '^## ' <slug>/lec.qmd` (+1 title) ≤ 30.
- Check no slide is overfull: rough rule ≤ ~10 lines of bullets or one figure
  plus ≤ 4 bullets, or code ≤ 12 lines plus ≤ 2 bullets.

## Report back (short)

Files written; slide count; proposed sidebar entries; glossary terms
(term — one-line definition); packages to add to the install page; anything
you were unsure about or deferred.
