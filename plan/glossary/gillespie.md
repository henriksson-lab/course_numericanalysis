**Copy number** — the whole number of molecules of a species in one cell.
**Propensity** — $a_j(\mathbf x)$: probability per unit time that reaction $j$ happens next in state $\mathbf x$.
**Stoichiometry vector** — the change in molecule counts caused by one firing of a reaction (a column of the stoichiometry matrix).
**Chemical master equation (CME)** — linear ODE system for the probabilities of all states of a reaction network; one equation per state.
**Birth–death process** — production at a constant rate and first-order decay; its steady state is Poisson.
**Gillespie algorithm (SSA, direct method)** — exact simulation of the CME: exponential waiting time with rate $a_0$, reaction chosen with probability $a_j/a_0$.
**Total propensity** — $a_0 = \sum_j a_j$, the rate at which any reaction happens.
**Telegraph (two-state) promoter model** — promoter switches between off and on; transcription only when on, giving bursts of mRNA.
**Transcriptional bursting** — mRNA made in bursts separated by silent periods, causing cell-to-cell variability.
**Burst size** — mean number of mRNA per burst, $k_m/k_\text{off}$ in the telegraph model.
**Fano factor** — variance divided by mean of a count; 1 for Poisson, > 1 for bursty production; depends on units.
**Constitutive expression** — transcription at a constant rate, giving Poisson mRNA numbers.
**Extinction (stochastic)** — reaching an absorbing state with zero individuals, which the ODE never reaches.
**Tau-leaping** — approximate SSA: in a fixed step $\tau$, each reaction fires a Poisson$(a_j\tau)$ number of times.
**Chemical Langevin equation** — SDE approximation of the CME: drift = ODE right-hand side, noise variance = sum of propensities.
**Euler–Maruyama method** — Euler's method for SDEs: an Euler step plus a normal kick with standard deviation proportional to $\sqrt h$.
**Stochastic differential equation (SDE)** — an ODE with a random noise term driven by Brownian motion $W$.

## Install

None (base R, numpy/scipy, and the Rust crates `rand`, `rand_distr` already listed).
