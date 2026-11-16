# Numerical analysis for biologists: course site sources

Teacher information. This file is not part of the published site.

Planning documents (teacher-only, not rendered):

- [audience.md](audience.md) — who the course is for and what we assume.
- [designprinciples.md](designprinciples.md) — how material is built.
- [contentplan.md](contentplan.md) — the planned units and their order.

## Layout

One folder per unit, each with `lec.qmd` (reveal.js deck), `exercise.qmd`
and `quiz.qmd`. The order of units lives in the sidebar in `_quarto.yml`, not
in folder names. Interactive demos are TypeScript in `demos/src/<name>/`
(see [demos/README.md](demos/README.md)). Figures made by a script live next to
the unit as `make-<unit>-figs.R` / `.py` and are regenerated with `make figs`.

## Building this site

```bash
npm i          # once: esbuild + typescript for the demos
make           # render the whole site to docs/ (then git add)
make render    # render only
make preview   # start a live-reload preview server
make clean     # delete docs/ and the Quarto cache
```

You will need [Quarto](https://quarto.org/docs/get-started/), a recent
[R](https://cran.r-project.org/), and Python 3 with numpy/scipy/matplotlib.
