// rust/examples/gradient_methods_c6.rs  —  cargo run --example gradient_methods_c6
// C6: argmin's Nelder-Mead, BFGS, L-BFGS and nonlinear CG on the 4PL; then the
// dose in nM, rescaled parameters, and log10(EC50). Like optim without `gr`,
// the gradient is by central differences with an absolute step 1e-3.
use argmin::core::{Executor, State, TerminationStatus};
use argmin::solver::conjugategradient::{beta::PolakRibiere, NonlinearConjugateGradient};
use argmin::solver::linesearch::MoreThuenteLineSearch;
use argmin::solver::neldermead::NelderMead;
use argmin::solver::quasinewton::{BFGS, LBFGS};
include!("gradient_methods_shared/defs.rs");

// p = q * scale; optionally p[2] = 10^p[2]; doses multiplied by `unit`
#[derive(Clone)]
struct Fd { unit: f64, scale: [f64; 4], log_ec50: bool }
impl Fd {
    fn par(&self, q: &[f64]) -> Vec<f64> {
        let mut p: Vec<f64> = q.iter().zip(&self.scale).map(|(q, s)| q * s).collect();
        if self.log_ec50 { p[2] = 10f64.powf(p[2]); }
        p
    }
}
impl CostFunction for Fd {
    type Param = Vec<f64>;
    type Output = f64;
    fn cost(&self, q: &Vec<f64>) -> Result<f64, Error> {
        let p = self.par(q);
        Ok(doses().iter().zip(Y.iter()).map(|(&c, &y)| (y - f4(&p, c * self.unit)).powi(2)).sum())
    }
}
impl Gradient for Fd {
    type Param = Vec<f64>;
    type Gradient = Vec<f64>;
    fn gradient(&self, q: &Vec<f64>) -> Result<Vec<f64>, Error> {
        (0..q.len()).map(|j| {
            let (mut a, mut b) = (q.clone(), q.clone());
            a[j] += 1e-3; b[j] -= 1e-3;
            Ok((self.cost(&a)? - self.cost(&b)?) / 2e-3)
        }).collect()
    }
}
fn eye() -> Vec<Vec<f64>> { (0..4).map(|i| (0..4).map(|j| (i == j) as i32 as f64).collect()).collect() }

fn run(name: &str, prob: Fd, q0: Vec<f64>, method: &str) {
    let ls = MoreThuenteLineSearch::new();
    let res = match method {
        "NM" => {
            let mut simplex = vec![q0.clone()];
            for j in 0..4 { let mut v = q0.clone(); v[j] += 0.1 * q0[j].abs().max(1.0); simplex.push(v); }
            let solver = NelderMead::new(simplex).with_sd_tolerance(1e-8).unwrap();
            Executor::new(prob.clone(), solver).configure(|s| s.max_iters(500)).run().map(|r| r.state().clone())
                .map(|s| (s.get_best_param().unwrap().clone(), s.get_best_cost(), s.get_iter(), s.get_func_counts().clone(), s.get_termination_status().clone()))
        }
        "BFGS" => Executor::new(prob.clone(), BFGS::new(ls))
            .configure(|s| s.param(q0.clone()).inv_hessian(eye()).max_iters(100)).run()
            .map(|r| { let s = r.state().clone(); (s.get_best_param().unwrap().clone(), s.get_best_cost(), s.get_iter(), s.get_func_counts().clone(), s.get_termination_status().clone()) }),
        "LBFGS" => Executor::new(prob.clone(), LBFGS::new(ls, 5))
            .configure(|s| s.param(q0.clone()).max_iters(100)).run()
            .map(|r| { let s = r.state().clone(); (s.get_best_param().unwrap().clone(), s.get_best_cost(), s.get_iter(), s.get_func_counts().clone(), s.get_termination_status().clone()) }),
        _ => Executor::new(prob.clone(), NonlinearConjugateGradient::new(ls, PolakRibiere::new()).restart_iters(10))
            .configure(|s| s.param(q0.clone()).max_iters(100)).run()
            .map(|r| { let s = r.state().clone(); (s.get_best_param().unwrap().clone(), s.get_best_cost(), s.get_iter(), s.get_func_counts().clone(), s.get_termination_status().clone()) }),
    };
    match res {
        Ok((q, cost, iters, counts, status)) => {
            let p = prob.par(&q);
            let why = match status { TerminationStatus::Terminated(r) => r.text().to_string(), _ => "running".into() };
            println!("{:<16} {:<5} SSE {:10.4}  EC50 {:9.4}  hill {:7.3}  iter {:3}  cost calls {:4}  ({})",
                     name, method, cost, p[2], p[3], iters, counts.get("cost_count").copied().unwrap_or(0), why);
        }
        Err(e) => println!("{:<16} {:<5} error: {}", name, method, e),
    }
}

fn main() {
    let um = Fd { unit: 1.0, scale: [1.0; 4], log_ec50: false };
    for m in ["NM", "BFGS", "LBFGS", "CG"] { run("µM", um.clone(), P0.to_vec(), m); }
    let nm = Fd { unit: 1000.0, scale: [1.0; 4], log_ec50: false };
    run("nM, EC50 = 1", nm.clone(), vec![0.0, 100.0, 1.0, 1.0], "BFGS");
    run("nM, EC50 = 1000", nm, vec![0.0, 100.0, 1000.0, 1.0], "BFGS");
    let sc = Fd { unit: 1000.0, scale: [1.0, 100.0, 1000.0, 1.0], log_ec50: false };
    run("nM, rescaled", sc, vec![0.0, 1.0, 1.0, 1.0], "BFGS");
    let lg = Fd { unit: 1000.0, scale: [1.0; 4], log_ec50: true };
    run("nM, log10 EC50", lg, vec![0.0, 100.0, 3.0, 1.0], "BFGS");
}
