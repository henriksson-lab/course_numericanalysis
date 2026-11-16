// C6: count function evaluations: Brent (argmin's BrentRoot), Newton, secant.
use argmin::core::{CostFunction, Error, Executor};
use argmin::solver::brent::BrentRoot;
use std::cell::Cell;
use std::sync::atomic::{AtomicU64, Ordering};

// global evaluation counter for the argmin run (argmin owns the function object)
static BRENT_EVALS: AtomicU64 = AtomicU64::new(0);

fn f(c: f64) -> f64 {
    0.3 * c / (0.1 + c) + 0.7 * c.powi(2) / (25.0 + c.powi(2)) - 0.5
}
fn df(c: f64) -> f64 {
    0.03 / (0.1 + c).powi(2) + 35.0 * c / (25.0 + c.powi(2)).powi(2)
}

fn newton(f: impl Fn(f64) -> f64, df: impl Fn(f64) -> f64, x0: f64, tol: f64, maxit: usize) -> f64 {
    let mut x = x0;
    for _ in 0..maxit {
        let step = f(x) / df(x);
        x -= step;
        if step.abs() <= tol * x.abs() {
            break;
        }
    }
    x
}

fn secant(f: impl Fn(f64) -> f64, mut x0: f64, mut x1: f64, tol: f64, maxit: usize) -> f64 {
    let (mut f0, mut f1) = (f(x0), f(x1));
    for _ in 0..maxit {
        if f1 == f0 {
            break;
        }
        let step = f1 * (x1 - x0) / (f1 - f0);
        x0 = x1;
        f0 = f1;
        x1 -= step;
        f1 = f(x1);
        if step.abs() <= tol * x1.abs() {
            break;
        }
    }
    x1
}

struct TwoSite;
impl CostFunction for TwoSite {
    type Param = f64;
    type Output = f64;
    fn cost(&self, c: &f64) -> Result<f64, Error> {
        BRENT_EVALS.fetch_add(1, Ordering::Relaxed);
        Ok(f(*c))
    }
}

fn main() {
    let c50 = 3.261305265055212;
    let report = |name: &str, evals: u64, root: f64| {
        println!("{:<22} evals {:3}   error {:.1e}", name, evals, (root - c50).abs());
    };
    for tol in [1e-4, 1e-12] {
        let res = Executor::new(TwoSite, BrentRoot::new(0.01, 100.0, tol))
            .configure(|s| s.max_iters(100))
            .run()
            .unwrap();
        let evals = BRENT_EVALS.swap(0, Ordering::Relaxed);
        report(&format!("BrentRoot, tol = {:.0e}", tol), evals, res.state.best_param.unwrap());
    }
    // a Cell is a counter that closures may change without `mut`
    let count = Cell::new(0u64);
    let fc = |c: f64| { count.set(count.get() + 1); f(c) };
    let dfc = |c: f64| { count.set(count.get() + 1); df(c) };
    let r = newton(fc, dfc, 1.0, 1e-12, 50);
    report("Newton from 1", count.replace(0), r);
    let r = secant(fc, 1.0, 2.0, 1e-12, 50);
    report("secant from 1, 2", count.replace(0), r);
}
