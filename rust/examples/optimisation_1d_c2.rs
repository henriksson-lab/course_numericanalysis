// C2: Brent's minimiser (argmin's BrentOpt) with its default and a tight tolerance,
// counting function evaluations.
use argmin::core::{CostFunction, Error, Executor, State};
use argmin::solver::brent::BrentOpt;
use std::sync::atomic::{AtomicU64, Ordering};

// evaluation counter (argmin owns the objective, so we count in a global)
static EVALS: AtomicU64 = AtomicU64::new(0);

const T: [f64; 7] = [0.0, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0];
const Y: [f64; 7] = [1.000, 0.804, 0.684, 0.484, 0.217, 0.094, 0.083];

struct Sse;
impl CostFunction for Sse {
    type Param = f64;
    type Output = f64;
    fn cost(&self, k: &f64) -> Result<f64, Error> {
        EVALS.fetch_add(1, Ordering::Relaxed);
        Ok(T.iter().zip(Y.iter()).map(|(t, y)| (y - (-k * t).exp()).powi(2)).sum())
    }
}

// minimise on [a, b]; tol = None keeps argmin's default (eps = sqrt(machine eps), t = 1e-5)
fn brent(a: f64, b: f64, tol: Option<f64>) -> Result<(f64, u64), Error> {
    EVALS.store(0, Ordering::Relaxed);
    let mut solver = BrentOpt::new(a, b);
    if let Some(t) = tol {
        solver = solver.set_tolerance(f64::EPSILON.sqrt(), t);
    }
    let res = Executor::new(Sse, solver).configure(|s| s.max_iters(100)).run()?;
    let st = res.state();
    let nev = EVALS.load(Ordering::Relaxed);
    Ok((*st.get_best_param().unwrap(), nev))
}

fn main() -> Result<(), Error> {
    let (k1, n1) = brent(0.0, 1.0, None)?;
    let (k2, n2) = brent(0.0, 1.0, Some(1e-10))?;
    println!("{:.10} {}", k1, n1);
    println!("{:.10} {}", k2, n2);
    println!("{:.6}", 2f64.ln() / k2); // half-life (h)
    let (k3, _) = brent(0.5, 2.0, Some(1e-10))?; // minimum outside the interval
    println!("{:.7}", k3);
    Ok(())
}
