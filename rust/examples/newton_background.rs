// Background: the two-site dose-response curve, a hand-written Newton, and a
// bracketing check with argmin's BrentRoot (Brent's method, as in uniroot/brentq).
use argmin::core::{CostFunction, Error, Executor};
use argmin::solver::brent::BrentRoot;

fn r(c: f64) -> f64 {
    0.3 * c / (0.1 + c) + 0.7 * c.powi(2) / (25.0 + c.powi(2))
}
fn f(c: f64) -> f64 {
    r(c) - 0.5
}
fn df(c: f64) -> f64 {
    0.03 / (0.1 + c).powi(2) + 35.0 * c / (25.0 + c.powi(2)).powi(2)
}

// returns (root, iterations); `impl Fn(f64) -> f64` accepts any function of one f64
fn newton(f: impl Fn(f64) -> f64, df: impl Fn(f64) -> f64, x0: f64, tol: f64, maxit: usize) -> (f64, usize) {
    let mut x = x0;
    for k in 1..=maxit {
        let step = f(x) / df(x);
        x -= step;
        if step.abs() <= tol * x.abs() {
            return (x, k);
        }
    }
    eprintln!("warning: no convergence after {maxit} steps");
    (x, maxit)
}

// argmin wants the function as a struct with a `cost` method
struct TwoSite;
impl CostFunction for TwoSite {
    type Param = f64;
    type Output = f64;
    fn cost(&self, c: &f64) -> Result<f64, Error> {
        Ok(f(*c))
    }
}

fn main() {
    let (root, it) = newton(f, df, 1.0, 1e-12, 50);
    println!("{}", root);
    println!("{}", it);
    let solver = BrentRoot::new(0.01, 100.0, 1e-12);
    let res = Executor::new(TwoSite, solver).configure(|s| s.max_iters(100)).run().unwrap();
    println!("{}", res.state.best_param.unwrap());
}
