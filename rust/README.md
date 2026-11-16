# Rust solutions

Every Rust snippet shown on an exercise or quiz page lives here as a runnable
example, so it is compiled and run, not just written.

```bash
cd rust
cargo run --example floating_point_c3     # one solution
cargo build --examples                    # check that all of them compile
```

File naming: `examples/<unit>_<task>.rs`, with the unit folder name in
snake_case and the task label in lower case (`intro_c4.rs`,
`error_conditioning_b2_hint.rs` for a hint).

Allowed crates (see `Cargo.toml`); use the standard library when it is enough:

| Crate | For | R / Python equivalent |
|--|--|--|
| `nalgebra` | dense linear algebra | `solve`, `%*%` / `numpy.linalg` |
| `sprs` | sparse matrices | `Matrix` / `scipy.sparse` |
| `argmin`, `argmin-math` | optimisation | `optim` / `scipy.optimize` |
| `statrs` | distributions, special functions | `dpois`, `lgamma` / `scipy.stats`, `scipy.special` |
| `rand`, `rand_distr` | seeded random numbers | `set.seed`, `rnorm` / `numpy.random` |
| `quadrature` | numerical integration | `integrate` / `scipy.integrate.quad` |
| `rustfft` | FFT | `fft` / `numpy.fft`, `scipy.fft` |
| `ode_solvers` | ODE solvers (Dopri5, Rk4, Dop853) | `deSolve::ode` / `scipy.integrate.solve_ivp` |
| `minilp` | linear programming | `lpSolve::lp` / `scipy.optimize.linprog` |
