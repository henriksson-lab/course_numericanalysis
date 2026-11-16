// C2: maximum-likelihood rate of a zero-truncated Poisson as the root of the score.
use argmin::core::{CostFunction, Error, Executor};
use argmin::solver::brent::BrentRoot;

const N: f64 = 100.0; // cells
const SX: f64 = 173.0; // total UMIs: 1*52 + 2*30 + 3*12 + 4*5 + 5*1

fn score(l: f64) -> f64 {
    SX / l - N / (-(-l).exp_m1()) // -expm1(-l) = 1 - exp(-l)
}
fn dscore(l: f64) -> f64 {
    -SX / l.powi(2) + N * (-l).exp() / (-l).exp_m1().powi(2)
}

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

struct Score;
impl CostFunction for Score {
    type Param = f64;
    type Output = f64;
    fn cost(&self, l: &f64) -> Result<f64, Error> {
        Ok(score(*l))
    }
}

fn main() {
    // (a) Newton from the mean, and the standard error
    let (lhat, it) = newton(score, dscore, SX / N, 1e-12, 50);
    println!("lambda-hat = {:.10} after {} steps, SE = {:.4}", lhat, it, 1.0 / (-dscore(lhat)).sqrt());
    // (b) many starting values
    for l0 in [0.2, 0.5, 1.0, 2.0, 2.5, 2.7, 2.8, 3.0, 4.0] {
        let mut l: f64 = l0;
        for _ in 0..50 {
            l -= score(l) / dscore(l);
        }
        println!("start {:.1} -> {:.6}", l0, l);
    }
    // (c) bracketing with Brent's method (argmin's BrentRoot)
    let res = Executor::new(Score, BrentRoot::new(0.01, 10.0, 1e-12))
        .configure(|s| s.max_iters(100))
        .run()
        .unwrap();
    println!("BrentRoot: {:.10} after {} iterations", res.state.best_param.unwrap(), res.state.iter);
}
