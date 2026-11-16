# Pass 2 brief: cross-references and consistency (teacher-only)

Pass 1 wrote every unit with cross-references as plain italic titles
(*Newton and secant methods*) and logged each one in `plan/xref/<slug>.md`.
Pass 2 turns them into links and makes sure that what a unit says about
another unit is true.

## Avoid permission prompts

Work only in the project folder and
`/big/henriksson/temp/claude-tmp/claude-1001/-data-henriksson-github-course-numericanalysis/f8db0def-e4f4-45f6-aad8-5433352f3d1d/scratchpad`.
No `/tmp` (use `mktemp -d -p <scratchpad>`), no git, no network, no edits
under `.claude/` or `.git/`, no `pkill`, no `rm -rf` outside the scratchpad.

## Read first

`plan/units.md` (titles), `designprinciples.md`, `plan/progress.md` (the
"Audit notes" list — items for your units are yours to fix), and for each of
your units `plan/xref/<slug>.md`.

## For each of your units

1. **Verify every cross-reference.** For each entry in the xref file (and any
   other mention of another unit you find in `lec.qmd`, `exercise.qmd`,
   `quiz.qmd`), open the target unit and check that the claim is true: the
   example, the numbers, the notation, the function names, what it covers.
   If it is not true, **fix the text in your unit** (never edit another
   worker's unit). Typical fixes: a recap slide that assumes a different
   example, a "as we saw in …" that the target doesn't actually show, a
   forward reference to content that ended up elsewhere.
2. **Turn titles into links.** First mention per page (and in "Where this
   leads") becomes a link: `[*Title*](../<slug>/lec.qmd)` — or
   `../<slug>/exercise.qmd#<anchor>` when pointing at a specific task.
   Links work on slides too. Later mentions on the same page/slide can stay
   plain italics. Do not link to units that don't exist; every slug is in
   `plan/units.md`.
3. **Running examples.** Where your unit reuses a running example defined in
   an earlier unit (see "Running examples defined here" in the xref files),
   make the numbers and names agree with the earlier unit. The earlier unit
   wins. If agreeing would require recomputing large parts of your unit and
   the difference would not confuse a student, instead add one plain sentence
   saying what differs and why (e.g. "here we use a cell-centred grid,
   because …"). Any code you change must be re-run in all three languages
   (R `Rscript --vanilla`, Python `python3 -I`, Rust
   `cd rust && cargo run -q --example <name>`) and the outputs on the page
   updated to the real ones; regenerate figures with the unit's
   `make-*-figs.R` if their data change.
4. **Audit-note items** in `plan/progress.md` that belong to your units and
   concern consistency or wording: fix them. Leave fact-checking against the
   literature (claims "from memory") for pass 3, but do not make them worse.
5. **Self-contained code tabs.** If a Python or R tab on an exercise page
   relies on functions defined only in another tab/file (not on the page),
   make it self-contained or add a one-line pointer to where the function is
   defined on the same page.
6. Update `plan/xref/<slug>.md`: mark each entry `[ok]`, `[fixed: …]` or
   `[note: …]`.

## Check

- `quarto render <slug>/lec.qmd`, `exercise.qmd`, `quiz.qmd` succeed (others
  render concurrently; retry on glob/site_libs errors).
- Deck stays ≤ 30 slides.
- `cd rust && cargo build -q --examples` has no warnings in your files.

## Report

Write `plan/pass2/<your group name>.md`: per unit, what you fixed, and any
problem that needs a change **in a unit outside your group** (say which unit,
which file/slide, what to change). Reply with ≤ 10 lines.
