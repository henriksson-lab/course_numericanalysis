**Pseudo-random number generator** — A deterministic algorithm that produces a long sequence of numbers that behave like independent uniform draws.
**Seed** — The starting state of a pseudo-random generator; the same seed (and software) gives the same stream.
**Linear congruential generator (LCG)** — The generator x_{n+1} = (a x_n + c) mod m; simple, short period, lattice structure.
**Period** — The number of steps after which a generator's sequence repeats.
**Mersenne Twister** — R's default generator (period 2^19937 − 1).
**PCG64** — numpy's default generator, used by `np.random.default_rng`.
**ChaCha12 / StdRng** — The default seeded generator of Rust's `rand` crate.
**Inverse-transform sampling** — Drawing X = F⁻¹(U) with U uniform, so that X has CDF F.
**Rejection sampling** — Drawing from a proposal g and accepting with probability f/(M g); accepted points follow f, the acceptance rate is 1/M.
**Monte Carlo integration** — Estimating an integral as (b − a) times the average of f at random points.
**Standard error (Monte Carlo)** — s/√N, the standard deviation of a Monte Carlo estimate; ±2 SE is an approximate 95% error bar.
**1/√N convergence** — The Monte Carlo error rate, independent of the dimension when the variance is finite.
**Infinite variance** — When ∫ f² diverges; the estimate still converges but slowly and in jumps, and the SE is unreliable.
**Curse of dimensionality** — The cost of a grid with n points per axis grows as n^d, so its error for a fixed budget N falls only like N^(−p/d).
**Marginal likelihood (evidence)** — Z = ∫ L(θ) π(θ) dθ, the normalising constant of the posterior; used in Bayes factors.
**Permutation-style (randomisation) test** — Computing a p-value by simulating the data many times under the null (e.g. shuffling genomic positions).
**Simulated p-value** — (k + 1)/(B + 1) for k of B simulations at least as extreme as the data.
**Importance sampling** — Drawing from a proposal g and weighting by w = f/g to estimate an integral with respect to f.
**Proposal distribution** — The distribution actually sampled from in rejection or importance sampling.
**Importance weight** — w(x) = f(x)/g(x), the correction for sampling from g instead of f.
**Effective sample size (ESS)** — (Σw)²/Σw², the number of equally weighted draws a weighted sample is worth.
**Rare event** — An event with probability so small that plain Monte Carlo rarely or never sees it.

## Install

- Nothing beyond the existing list. R: base and stats only. Python: numpy, scipy. Rust: `rand`, `rand_distr`, `statrs` (already in Cargo.toml).
