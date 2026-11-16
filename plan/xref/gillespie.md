# Cross-references from `gillespie` (Stochastic simulation)

- gillespie/lec.qmd:34 (notes, About this lecture) -> monte-carlo (SSA as a Monte Carlo method; standard errors) [ok]
- gillespie/lec.qmd:39 (slide "Recap: the mRNA ODE") -> euler (mRNA model m' = alpha - delta m from slide "ODEs as models") [fixed: linked; slide now states euler writes m' = alpha - delta m with alpha = 10/min, delta = 0.2/min (steady state 50), here k = 10/h, delta = 1/h]
- gillespie/lec.qmd:53 (notes, Recap) -> euler (there alpha = 10 per min, delta = 0.2 per min, steady state 50) [fixed: notes explain the renaming alpha -> k and the change of units]
- gillespie/lec.qmd:109 (notes, "The chemical master equation, in words") -> power-iteration (stationary distribution of a Markov chain as eigenvector for eigenvalue 1/0) [fixed: power-iteration uses eigenvalue 1 of a transition matrix; notes now say so (eigenvalue 0 of the rate matrix)]
- gillespie/lec.qmd:120 (notes, "Gillespie's idea") -> monte-carlo (inverse-transform sampling, exponential and discrete/base-composition example) [ok]
- gillespie/lec.qmd:211 (notes, "Many cells") -> monte-carlo (Monte Carlo standard error sd/sqrt(N); streams differ between languages) [ok]
- gillespie/lec.qmd:222 (notes, "Runs, mean and distribution") -> monte-carlo (testing on a known answer: binomial overlap count) [ok]
- gillespie/lec.qmd:336 (slide "When the mean is not the ODE: extinction") -> euler (SIR running example N = 10 000, beta = 0.3, gamma = 0.1) [fixed: linked; N, beta, gamma match euler; I0 = 1 here vs 10 there is stated]
- gillespie/lec.qmd:356 (notes, extinction) -> phase-plane (noise switching the toggle switch between stable states) [ok]
- gillespie/lec.qmd:367 (notes, "Single outbreaks vs the ODE") -> nonlinear-least-squares (fitting an ODE to a stochastic outbreak; "Part 9") [fixed: Part 9 fits deterministic models only; wording no longer implies NLS treats outbreaks]
- gillespie/lec.qmd:408 (notes, tau-leaping) -> euler (tau-leaping is to the SSA what Euler is to the ODE) [ok]
- gillespie/lec.qmd:433 (notes, Euler–Maruyama) -> euler (order 1 of Euler for ODEs) [ok]
- gillespie/lec.qmd:458 (Summary, "Next") -> optimisation-1d (start of Part 9) [fixed: linked]
- gillespie/lec.qmd:461 (notes, Summary) -> nonlinear-least-squares, monte-carlo [ok]
- gillespie/exercise.qmd:135 (B7(c)) -> euler (SIR example starts with I = 10) [fixed: linked; I = 10 matches euler]
- gillespie/exercise.qmd:921 (Solution C6) -> monte-carlo (rare-event problem) [fixed: linked]
- gillespie/exercise.qmd:926 (Where this leads) -> optimisation-1d, nonlinear-least-squares (fitting rates / ODE parameters) [fixed: linked]
- gillespie/exercise.qmd:927 (Where this leads) -> identifiability (burst size and frequency identifiable, not all four telegraph rates) [fixed: linked; note: identifiability does not treat the telegraph model, sentence is a forward pointer to the topic]
- gillespie/exercise.qmd:928 (Where this leads) -> reaction-diffusion (stochastic reaction–diffusion with the SSA) [fixed: linked; added that reaction-diffusion treats deterministic PDEs only]

## Running examples defined here

- **Low-expressed gene, birth–death** (time in hours): k = 10 per h, delta = 1 per h, m(0) = 0,
  ODE m(t) = 10(1 - e^{-t}); steady state Poisson(10). R set.seed(1), 2000 cells with `ssa()` at
  t = 0.5, 1, 2, 10: means 3.8575, 6.3285, 8.664, 10.1185 (SE 0.044–0.070), var at 10 h 9.749.
  Per-event average bug: 10.5 (exact) vs time average 10.
- **Constitutive gene with protein**: k = 10, delta_m = 1, k_p = 2 per mRNA, delta_p = 0.2 (per h);
  means 10 and 100; exact Fano factors 1 (mRNA) and 1 + k_p/(delta_m + delta_p) = 2.667 (protein).
- **Telegraph (bursty) gene**: k_on = 0.5, k_off = 4.5, k_m = 100 (per h), same delta_m, k_p, delta_p;
  fraction on 0.1, mean mRNA 10, protein 100 (same ODE as constitutive); burst size k_m/k_off = 22.2;
  exact Fano mRNA 16, protein 32.47 (moment/Lyapunov equations). R set.seed(1), 1000 cells at t = 30 h:
  mean_m 10.107, Fano_m 15.26, mean_p 100.21, Fano_p 32.93; 18.9% of cells with zero mRNA.
- **SIR extinction**: euler's SIR (N = 10 000, beta = 0.3, gamma = 0.1) from I0 = 1: P(extinct) = 1/R0 = 1/3
  (R set.seed(1), 10 000 runs stopping at I = 0 or I = 100: 0.3348, SE 0.0047); I0 = 2: 1/9, I0 = 3: 1/27.
- **Tau-leap variance** for birth–death: 2k/(delta(2 - delta tau)) = 20/(2 - tau).
