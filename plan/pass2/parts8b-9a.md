# Pass 2 report: parts8b-9a (stiffness, phase-plane, gillespie, optimisation-1d)

All italic cross-references to other units are now links (first mention per slide body in the decks, first
mention per page on exercise/quiz pages, and every item in "Where this leads"). The notes stay plain text.
Each entry in `plan/xref/<slug>.md` is marked `[ok]` / `[fixed: …]` / `[note: …]`.
No code outputs changed, so nothing was re-run and no figures were regenerated. No Rust files were touched.

## stiffness

- Overflow on "Try it: the fast eigenvalue sets the step" (741/720 px) is fixed. The demo `div` now has
  `style="max-width: 1120px; margin: 0 auto;"`; the SVG plots scale with the width, so the content is the same, only a
  little smaller. A headless re-probe gives 720 px. No other slide in the deck is over 720.
- Checked against the source units:
  - euler: drug clearance −0.8 C, factors 0.6/−0.6/−1.4, −14/19.6; 24-h oscillator, Euler 2850.
  - newton-systems: `newton_sys`, finite-difference Jacobian.
  - runge-kutta: ode45 = Dormand–Prince.
  - finite-differences: the one-sided second-order formula.
  - lu: 2/3 n³.
  - diffusion-implicit: stiff method of lines, a tridiagonal solve per step.
  
  All are correct.
- 12 titles are now links (lec 5, exercise 7).

## phase-plane

- The toggle switch numbers match newton-systems exactly: α = 3, n = 2, states (2.618034, 0.381966) /
  (1.213412, 1.213412), Newton basins 21/57/21 %, det J = −0.4186 at the middle state, exercise C3 23/54/23 %, and
  the singular J at α = 2 (39 steps).
- Repressilator: Sel'kov stays. The slide "Oscillations: glycolysis in yeast" now names the repressilator of
  *Runge–Kutta and adaptive step size* (linked, six variables) as a larger model that is analysed the same way. The
  notes give α = 216 and period ≈ 8.7, and say that Sel'kov is used because two variables fit in one phase plane.
  The exercise Background has a matching sentence with a link.
- Exercise E3 said that newton-systems shows Newton basins for the *biased* switch (21.8/57/21.2 %). It does not;
  newton-systems only has the symmetric switch. E3 now points to this lecture's slide for the biased numbers and
  quotes 21/57/21 % for the symmetric switch in newton-systems. The demo intro now says "the demo in *Newton for
  systems…* (symmetric switch only)".
- Self-contained tabs: the C2, C3, C4 and C6 R/Python tabs now have a one-line comment naming where
  `toggle`/`J`/`which_state`/`fd_jacobian`/`states2` are defined. A "solutions build on each other" paragraph
  already existed.
- Links: lec 5 + 1 new, exercise 6 + 2.

## gillespie

- The gene parameters are now stated explicitly on the recap slide. In euler the model is m' = α − δm with α = 10/min
  and δ = 0.2/min (steady state 50). Here the transcription rate is called k, with k = 10/h and δ = 1/h (steady
  state 10). The notes explain the renaming and the change of units.
- SIR checked: N = 10 000, β = 0.3, γ = 0.1 and R0 = 3 match euler; euler starts from I = 10, gillespie from I = 1.
  Both are stated. Exercise B7(c) "I = 10 (start of the SIR example)" is correct.
- Notes on the CME: power-iteration uses eigenvalue **1** of a discrete-time transition matrix. The text said "eigenvalue
  0 as for the Markov chains in Power iteration", which is now made precise.
- "Fitting an ODE to a stochastic outbreak … (Part 9)": Part 9 does not treat outbreaks. The wording now says only
  that the fitting methods are in Part 9.
- In "Where this leads", the reaction-diffusion item now says that the unit uses deterministic PDEs (it has no SSA).
  The identifiability item (burst size and frequency versus four telegraph rates) stays as a forward pointer to
  the topic. Identifiability itself does not discuss the telegraph model.
- Links: lec 2 + recap, exercise 6.

## optimisation-1d

- The zero-truncated Poisson numbers match newton exactly: counts 52/30/12/5/1, Σ = 173, λ̂ = 1.2184695258, Newton
  converges in 6 steps from the mean and fails from 3, `uniroot` on [0.01, 10]. The SE from nll'' = 56.9 is 0.1326,
  against newton's 0.1325.
- Other references checked and correct: the `<<-` counter and Rust `AtomicU64` from newton exercise C6; the logistic
  OD data and the 3-parameter MLE r = 0.7954 (finite-differences); line search (gradient-methods); periodogram
  (fourier-fft).
- Fixed: the lec notes (Summary) and the exercise "Where this leads" said that the mRNA decay fit "returns with two
  parameters" in *Nonlinear least squares and ODE fitting*. NLS fits a three-parameter PK curve and ODE models; it
  has no decay fit. The text is reworded.
- Links: lec 3, exercise 8, quiz 1.

## Renders

Every edited page rendered without errors: all four decks, all four exercises and the optimisation-1d quiz. Glob and
site_libs errors were retried. The stiffness, phase-plane and gillespie quizzes were not changed. A headless overflow
probe of all four decks finds no slide over 720 px. Slide counts: stiffness 26, phase-plane 26, gillespie 27,
optimisation-1d 28. The coordinator's extra item (optimisation-1d lec notes about line 482, the false "decay fit
returns with two parameters in NLS") is done; see the optimisation-1d section.

## Needed outside this group

- **nonlinear-least-squares** (optional): `plan/xref/optimisation-1d.md` had proposed reusing the mRNA decay series
  (t = 0…8 h, k̂ = 0.37539) as a two-parameter example. NLS did not take it up. No change is needed; this unit
  no longer promises it.
- **newton-systems** (optional, for the phase-plane comparison): basins exist only for the symmetric switch.
  phase-plane is now consistent with that, so no change is needed.
- Housekeeping, already logged by the pass-1 stiffness worker: a stray Chrome profile directory at
  `/big/henriksson/temp/claude-tmp/claude-1001/chr1`, outside the scratchpad. It was not touched.
