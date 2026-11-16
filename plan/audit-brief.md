# Pass 3 brief: correctness audit (teacher-only)

The units are written (pass 1) and cross-linked (pass 2). Rust outputs have
already been checked automatically (all match), and no slide overflows. This
pass checks everything else that can be wrong.

## Avoid permission prompts

Work only in the project folder and
`/big/henriksson/temp/claude-tmp/claude-1001/-data-henriksson-github-course-numericanalysis/f8db0def-e4f4-45f6-aad8-5433352f3d1d/scratchpad`.
No `/tmp` (use `mktemp -d -p <scratchpad>`), no git, no network, no edits
under `.claude/` or `.git/`, no `pkill`, no `rm -rf` outside the scratchpad.
Do not fact-check literature claims (that needs the web and is listed in
`plan/factcheck.md` for later); just don't make them worse.

## For each of your units

1. **R and Python code.** Extract every R and every Python code block from
   `lec.qmd` (R only), `exercise.qmd` and `quiz.qmd`, in page order. Blocks on
   one page may build on earlier blocks of the same language: run them
   cumulatively per page and language (one R session / one Python process per
   page), with `Rscript --vanilla` and `python3 -I`. Compare every printed
   result with what the page shows (output comments, `# Output:` blocks or
   output fences). For each mismatch: if the code is right, update the shown
   output to the real one; if the code is wrong, fix the code (and the Rust
   twin in `rust/examples/` if it shares the bug — re-run it). Machine-
   dependent timings may differ: accept them if the page says they are
   machine-dependent; otherwise add that remark. Random numbers: outputs must
   match for the stated seed.
2. **Figures.** Run `make-<slug>-figs.R` from the unit folder. It must run
   without errors, and every figure referenced in `lec.qmd`/`exercise.qmd`
   must exist. If a figure's data changed in pass 2 (see `plan/pass2/*.md`),
   check that the figure still agrees with the text.
3. **By-hand work and quiz answers.** Re-do the arithmetic of every "By hand"
   solution and every quiz answer. Check the reasoning, not just the numbers.
   Fix errors.
4. **Text versus numbers.** Where slide or exercise text quotes a number that
   comes from code (e.g. "converges in 6 steps", "error 3.2e-9"), check it
   against the real output.
5. **Items for your units** in `plan/progress.md` "Audit notes" that are not
   literature fact-checks: fix them or explain in your report why not.

## Check

- All three pages of each unit render (`quarto render <slug>/<page>.qmd`;
  retry on glob/site_libs errors — others render concurrently).
- Decks stay ≤ 30 slides; do not add slides.
- If you touched Rust: `cd rust && cargo build -q --examples` has no warnings
  in your files and the shown output is the real one (debug build, unless the
  header says `--release`).

## Report

Write `plan/audit/<your group name>.md`: per unit, what you checked and what
you fixed (one line each), anything you could not resolve. Reply with ≤ 10
lines.
