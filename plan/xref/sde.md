# Cross-references from `sde` (Stochastic differential equations)

- sde/lec.qmd:34 (notes, About this lecture) -> gillespie (one Euler–Maruyama slide; SDE theory listed as not covered there) [ok]
- sde/lec.qmd:39 (slide "Recap: noise beyond molecule counts") -> gillespie (chemical Langevin equation dm = (k - delta m) dt + sqrt(k + delta m) dW, EM one-liner from slide "Euler–Maruyama for stochastic differential equations") [ok]
- sde/lec.qmd:56 (notes, Recap) -> gillespie (states EM strong order 1/2, weak order 1; exercise "Not covered" line lists SDE theory) [ok]
- sde/lec.qmd:90 (notes, "Why steps of size sqrt h") -> diffusion-explicit (diffusion length sqrt(2 D t), slide "The diffusion time scale L^2/D") [ok]
- sde/lec.qmd:125 (notes, "D from the MSD in R") -> least-squares-qr (lm as least squares; here one column through the origin) [ok]
- sde/lec.qmd:150 (slide "The SDE and the Euler–Maruyama step") -> euler (Euler's method; g = 0 gives Euler) [ok]
- sde/lec.qmd:196 (notes, OU test problem) -> stiffness (test equation y' = lambda y) [ok]
- sde/lec.qmd:222 (notes, "EM against the exact recursion") -> stiffness (Euler stability limit |1 - theta h| < 1) [ok]
- sde/lec.qmd:239 (slide "Strong and weak error") -> convergence-cost (order read off as slope on a log-log plot) [ok]
- sde/lec.qmd:394 (notes, Failure 1) -> gillespie (8% NaN with k = 0.5 in the gillespie slide; here 6.4%, other random stream) [ok]
- sde/lec.qmd:397 (slide "Failure 2") + :419 (notes) -> euler, stiffness (drug clearance C' = -0.8 C, Euler stable for h < 2.5 h: stiffness slide "Recap: Euler's stability limit") [ok]
- sde/lec.qmd:424 (slide "The fix: drift-implicit Euler–Maruyama") -> stiffness (backward Euler) [ok]
- sde/lec.qmd:440 (slide "From paths to a PDE: Fokker–Planck") -> advection-diffusion (flow-plus-diffusion equation, velocity u and D) [ok]
- sde/lec.qmd:478 (notes, "Fokker–Planck in R") -> advection-diffusion (C + 2r <= 1, cell Péclet < 2 for central, upwind numerical diffusion u dx/2; "refine until P well below 2 and use central") [ok]
- sde/lec.qmd:512 (notes, Cost) -> monte-carlo (1/sqrt(N) law, slide "The error falls like 1/sqrt(N)") [ok]
- sde/lec.qmd:538 (Summary, "Next") -> optimisation-1d (start of Part 9) [ok]
- sde/exercise.qmd:204 (Hint 6) -> advection-diffusion (upwind numerical diffusion |u| dx/2) [ok]
- sde/exercise.qmd:1153 (Where this leads) -> optimisation-1d, nonlinear-least-squares (fitting model parameters; Part 9) [ok]
- sde/exercise.qmd:1154 (Where this leads) -> monte-carlo (Monte Carlo error; importance sampling for rare events) [ok]
- sde/exercise.qmd:1155 (Where this leads) -> advection-diffusion (scheme used for Fokker–Planck) [ok]
- sde/quiz.qmd:154 (Q12) -> phase-plane (toggle switch: two stable states with basins), monte-carlo (rare events) [ok]

## Running examples defined here

- **Bead in water / microchannel**: radius 0.5 µm, 20 °C, eta = 1.0 mPa s, D = 0.43 µm²/s = 430 nm²/ms; 10 frames/s, 100-frame tracks.
- **Optical trap (OU)**: x in nm, t in ms; theta = 1 per ms, D = 430 nm²/ms, sigma = sqrt(2D) = 29.3 nm/sqrt(ms); stationary variance 430 nm² (sd 20.7 nm). Release experiment: X(0) ~ N(60, 5²) nm; exact at 1 ms mean 22.07, variance 375.2.
- **GBM numerics test**: r = 0.8 per hour (OD running example rate), sigma = 0.5, N0 = 1, T = 1; R set.seed(1), M = 10 000, 2^10 fine steps: EM strong 0.0877/0.0400/0.0199/0.0099 at h = 1/16, 1/64, 1/256, 1/1024; Milstein 0.0565/0.0146/0.0037/0.0009.
- **GBM population (Itô)**: r = 0.1 per year, sigma = 0.6, N0 = 1000; median e^{-0.08 t}; 50 years: mean 148 413, median 18.3, P(N < 10) = 0.443.
- **Noisy drug clearance**: dC = -0.8 C dt + C dW (per hour), C0 = 10; EM mean-square stable only for h < 0.9375 h (ODE: h < 2.5 h).
