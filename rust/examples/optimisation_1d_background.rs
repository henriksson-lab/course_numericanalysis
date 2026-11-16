// Background: minimise the decay-rate SSE with Brent's method (argmin's BrentOpt,
// Rust's counterpart of R's optimize and scipy's minimize_scalar).
use argmin::core::{CostFunction, Error, Executor, State};
use argmin::solver::brent::BrentOpt;

const T: [f64; 7] = [0.0, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0]; // hours
const Y: [f64; 7] = [1.000, 0.804, 0.684, 0.484, 0.217, 0.094, 0.083]; // fraction left

// argmin wants the objective as a type with a `cost` method
struct Sse;
impl CostFunction for Sse {
    type Param = f64;
    type Output = f64;
    fn cost(&self, k: &f64) -> Result<f64, Error> {
        Ok(T.iter().zip(Y.iter()).map(|(t, y)| (y - (-k * t).exp()).powi(2)).sum())
    }
}

fn main() -> Result<(), Error> {
    // BrentOpt stops when the bracket is below eps*|x| + t: here eps = sqrt(machine eps), t = 1e-10
    let solver = BrentOpt::new(0.0, 1.0).set_tolerance(f64::EPSILON.sqrt(), 1e-10);
    let res = Executor::new(Sse, solver).configure(|s| s.max_iters(100)).run()?;
    let k = *res.state().get_best_param().unwrap();
    println!("{:.7}", k);                          // decay rate (per hour)
    println!("{:.9}", res.state().get_best_cost()); // SSE at the minimum
    println!("{:.6}", 2f64.ln() / k);               // half-life (hours)
    Ok(())
}
