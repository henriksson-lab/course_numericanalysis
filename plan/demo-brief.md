# Brief for writing an interactive demo (teacher-only, not rendered)

Read first: `demos/README.md` (rules), `demos/src/lib/ui.ts` (helpers: slider,
button, checkbox, svgEl, rng, auto-mount), **our own demos as the primary
reference** — `demos/src/stepsize/main.ts` with its Explore section in
`convergence-cost/exercise.qmd` and its slide in `convergence-cost/lec.qmd`,
and `demos/src/taylor/` — and, secondarily, one demo from the
sister course (read-only):
`/data/henriksson/github/course/quantbio-3MB00A/demos/src/pseudocount/main.ts`
with its embedding in `.../scrna-log/exercise.qmd` (search `data-demo`) and
`.../scrna-log/lec.qmd`. Also skim `.../demos/src/compensation/main.ts`.

## What to produce

1. `demos/src/<name>/main.ts` — plain TypeScript, DOM + SVG, no libraries, no
   network. Auto-mount via the helper in `lib/ui.ts`. Fixed seeds. Fits in
   about 1150×560 px; labels ≥ 18 px; preset/reset buttons; nothing hover-only;
   never bind arrow keys or space. Read settings from `data-*` attributes.
   Use the shared `demos/demo.css` classes; put demo-specific styles inline or
   in the TS. Colour-blind-safe colours (Okabe–Ito), consistent with the
   unit's figures.
2. **Exercise page**: add an `## Explore {#explore}` section to
   `<unit>/exercise.qmd`, placed after "In code" and before "Where this leads"
   (or wherever the reference puts it), with the embed block

   ````md
   ```{=html}
   <link rel="stylesheet" href="../demos/dist/demo.css">
   <div data-demo="<name>" ...></div>
   <script src="../demos/dist/<name>.js"></script>
   ```
   ````

   followed by 3–5 numbered tasks (prose only — no code tabsets in Explore) **E1, E2, …** that make the student *use*
   the demo to discover the unit's key idea, each with a collapsible
   Answer callout (same style as the page's other solutions). Answers quote
   numbers you actually read off the demo.
3. **Lecture**: add **one** slide to `<unit>/lec.qmd` that embeds the same
   demo (with `::: notes` describing what to show in the recording), at the
   point in the deck where it fits best. The deck must stay ≤ 30 slides
   including the title. Keep the static figure that the demo illustrates.

Numbers: JavaScript doubles are IEEE doubles, so floating-point behaviour is
the same as in R/numpy, but `Math.exp`/`Math.log` may differ in the last bit;
do not promise bit-identical agreement with R output on the page.

## Check before reporting

- `demos/build.sh <name>` builds without warnings; `make typecheck` passes
  (if it fails only in another demo's folder, say so and ignore it).
- `quarto render <unit>/exercise.qmd` and `quarto render <unit>/lec.qmd`
  succeed.
- Load the rendered exercise page in headless Chrome
  (`google-chrome --headless=new --no-sandbox --window-size=1280,900
  --dump-dom file://.../docs/<unit>/exercise.html`)
  and confirm the demo mounted (its SVG is present in the DOM, no visible
  error). Use a fresh `--user-data-dir=$(mktemp -d)` per Chrome call and do
  NOT use `--virtual-time-budget` (it hangs here); wrap in `timeout 60`.
  Screenshots via `--screenshot` are optional; retry once on an empty result.
- Other workers render and build concurrently; retry a failed `quarto render`
  once before investigating.
- Only touch: `demos/src/<name>/`, `demos/dist/<name>.js` (built),
  `<unit>/exercise.qmd`, `<unit>/lec.qmd`. Do not edit `_quarto.yml`, the
  shared lib, or other units — report needed changes instead.

## Report back (short)

Demo name and what it shows (also give a short sidebar label, 2–5 words, e.g. "choose the step size"); where it was embedded (slide title, exercise
section); slide count after the change; the E-tasks in one line each;
anything unsure.
