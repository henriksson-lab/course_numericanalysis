# Pass 3 audit: part10a (units 31–34)

Group: diffusion-explicit, diffusion-implicit, steady-state-2d, advection-diffusion.

How I ran the code:
- Every R and Python block of lec, exercise and quiz was run cumulatively, with one `Rscript --vanilla` / `python3 -I` process per page and language.
- Each printed line was compared with the output shown on the page, in both directions.
- Quiz pages have no code.
- Some lec slides use helpers that are defined elsewhere (`thomas`, `exact`, `embryo_steady`, `Tout`, `cobs`). I re-ran those slides with the definitions taken from the exercise Background or the figure script.
- All four `make-<slug>-figs.R` scripts run without errors. Every referenced figure exists.
- No figure data changed in pass 2.
- No Rust was touched.

## diffusion-explicit (31)

- Exercise code: all 12 R/Python blocks now match.
- **Fixed:** the C2 solution looped `for (D in …) for (L in …)` (Python: `for D …`, `for L …`), which overwrote the Background's `D` and `L`. Anyone running C3 and C4 after C2 in the same session got wrong numbers (dx = 50 µm, dt = 10 s, t D/L² = 0.033). The loop variables are now `Dj`/`Lj` in R and Python. The output is unchanged and was re-verified. Rust is not affected.
- Lec code: all 4 R blocks match.
- By hand B1–B8 and quiz Q1–Q11 re-done, all correct. This includes the B4 tables, B5 (480 steps), B6 eigenvalues and rates, and B8 counts.
- Text vs numbers checked, all correct:
  - slope −0.0030/µm at x = 100 µm, t = 30 s
  - 1611 steps / 64.4 s
  - 1 % criterion takes 1.6× longer
  - zigzag invisible for about 30 steps and clear by 60
  - rates 0.0395 to 40 per s, τ = 25.3 s
- progress.md note "slides not visually checked": no action needed. The brief states that no slide overflows.

## diffusion-implicit (32)

- Exercise code: all 14 R/Python blocks match.
- Lec code: all R blocks match (with the Background helpers):
  - 15 BE steps: 0.8781…
  - order ratios 1.92/4.03…
  - T200 sequence and −9.05 3.98 4.00
  - lsoda 128 / bdf 85 / ode45 3563 steps
- By hand B1–B6 re-done (Thomas fractions, growth factors, cost counts, embryo τ/λ, PCR wall t₁/Bi), all correct. The exact chamber values (0.480, 0.157, 0.034) and (0.617, 0.315, 0.121) were confirmed with the series.
- Quiz Q1–Q11 all correct. Q4's 1.6 / 0.36 / 1.46 was reproduced, and r = 4 does stay in [0, 1].
- progress.md notes:
  - Self-contained tabs: settled. Pass 2's "uses … from the Background" comments are the convention, and the cumulative run confirms every tab works after the Background.
  - h = 10 W/(m²K): already hedged in the text (range 5–25).
  - [M1], [M2], [H1], [H2]: literature, left for factcheck.md.

## steady-state-2d (33)

- Exercise code: all 14 blocks match. Lec code: all 6 blocks match.
- Other lec numbers re-computed with the figure script's model, all correct:
  - the fill-in table (240 317 / 1 984 637 / 16 129 277 natural; 43 825 / 257 722 / 1 475 954 reordered)
  - 9.5 million dense entries
  - plain Newton 12 steps to −27.61
  - clipped 8 steps to 0.564
  - start from 0: 12 steps, same field to 4e-9
  - constant uptake −27.71 / 86.15
  - 5.7 / 10.0 / 16.3 % at 5 / 10 / 20 µM
  - E2 table (32.2, 12.4, 8.24)
- **Threshold density (fig10 notes), verified with the model:** Newton with Michaelis–Menten, h = 50 µm, root-found on the density.
  - The minimum first reaches 20 / 10 / 5 µM at **1.62 / 1.73 / 1.80 e5** cells/cm².
  - The hand estimate "5 µM at about 1.7e5" was low.
  - The notes now give the model values. They also give the constant-uptake estimate (1.58 / 1.67 / 1.71 e5) and say why it comes out slightly early.
- **Fixed:** in the B1(a) solution, the derivation said "add to twice the first: c₂ = 337.5 − 200". The result was right but the step was wrong. It now reads: adding the first equation gives c₂ = 168.75 − 31.25.
- The other by-hand items (B2–B8: trench, numbering, Kronecker, Newton iterates on both roots, fill-in estimates) and quiz Q1–Q11 are all correct.
- progress.md notes: [B1] "mechanics not nutrients" is literature (factcheck). Grid conventions were done in pass 2.

## advection-diffusion (34)

- Exercise code: all 14 blocks match.
- Lec code: all R blocks match, including march 0.3944/0.1361 and optimize 398.1985. Brent needs 19 evaluations, as the text says.
- Figure-script printouts agree with the text: 0.4 % error at dx = 50 µm, zigzag 8e39, fractions 47 / 19 / 0.1 % at 2 mm, staircase k/7.
- By hand B1–B6 re-done, all correct: Pe and mixing lengths, the N = 4 central and upwind solutions and ρ, D_num, the CFL weights, the 168 000 steps, the Christmas tree and B6.
- **Fixed:** B6(c) said "398 · 0.8 = 318.6". It now reads 398.2 · 0.8.
- **Fixed:** C6 answer (b) said the fit "settles near 396 only for Δx ≤ 25 µm", but Δx = 25 µm gives 398.2. It now says that from Δx ≤ 25 µm the estimate is within 1 % of the converged value, about 396.
- Quiz Q1–Q11 all correct (Q7 168 000 steps, Q8 1.6 mm / 17 %, Q10 505 / 396 / 27 %).
- progress.md notes:
  - Grid, mixing length and the E4 stability condition were done in pass 2.
  - Code blocks of 13–14 lines are a style point. I did not split them (that would risk overflow).
  - Kamholz & Yager is left for factcheck.

## Checks

- Pages render:
  - diffusion-explicit: lec, exercise, quiz
  - diffusion-implicit: lec, exercise, quiz
  - steady-state-2d: lec, exercise, quiz
  - advection-diffusion: lec, exercise, quiz

  All 12 pages rendered. Some needed 2–4 tries because of concurrent glob/site_libs errors.
- Decks: 25 / 25 / 27 / 28 `##` slides. No slides added.

## Unresolved

- None.
