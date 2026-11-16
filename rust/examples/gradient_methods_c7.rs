// rust/examples/gradient_methods_c7.rs  —  cargo run --example gradient_methods_c7
// C7: the growth curve of "Why numerical analysis?", time in minutes vs hours,
// with argmin's Nelder-Mead and BFGS (central-difference gradient, step 1e-6).
use argmin::core::{CostFunction, Error, Executor, Gradient, State};
use argmin::solver::linesearch::MoreThuenteLineSearch;
use argmin::solver::neldermead::NelderMead;
use argmin::solver::quasinewton::BFGS;

const OD: [f64; 13] = [0.001, 0.047, 0.139, 0.282, 0.581, 0.859, 1.049,
                       1.161, 1.152, 1.217, 1.182, 1.176, 1.185];
fn logistic(p: &[f64], t: f64) -> f64 {        // p = (r, K, N0)
    p[1] / (1.0 + (p[1] / p[2] - 1.0) * (-p[0] * t).exp())
}
#[derive(Clone)]
struct Growth { t: Vec<f64> }
impl CostFunction for Growth {
    type Param = Vec<f64>;
    type Output = f64;
    fn cost(&self, p: &Vec<f64>) -> Result<f64, Error> {
        Ok(self.t.iter().zip(OD.iter()).map(|(&t, &y)| (y - logistic(p, t)).powi(2)).sum())
    }
}
impl Gradient for Growth {
    type Param = Vec<f64>;
    type Gradient = Vec<f64>;
    fn gradient(&self, p: &Vec<f64>) -> Result<Vec<f64>, Error> {
        (0..3).map(|j| {
            let h = 1e-6 * p[j].abs().max(1.0);
            let (mut a, mut b) = (p.clone(), p.clone());
            a[j] += h; b[j] -= h;
            Ok((self.cost(&a)? - self.cost(&b)?) / (2.0 * h))
        }).collect()
    }
}

fn main() -> Result<(), Error> {
    let t_min: Vec<f64> = (0..=12).map(|i| 120.0 * i as f64).collect();
    let start = vec![0.5, 1.0, 0.05];
    // (a)
    println!("{:.6?}", Growth { t: t_min.clone() }.gradient(&start)?);
    // (b)
    for (unit, t) in [("minutes", t_min.clone()), ("hours", t_min.iter().map(|t| t / 60.0).collect())] {
        let prob = Growth { t };
        let mut simplex = vec![start.clone()];
        for j in 0..3 { let mut v = start.clone(); v[j] += 0.1; simplex.push(v); }
        let r = Executor::new(prob.clone(), NelderMead::new(simplex).with_sd_tolerance(1e-10)?)
            .configure(|s| s.max_iters(500)).run()?;
        let p = r.state().get_best_param().unwrap();
        println!("{:<7} Nelder-Mead r {:.6}  K {:.4}  SSE {:.6}", unit, p[0], p[1], r.state().get_best_cost());
        let eye = vec![vec![1.0, 0.0, 0.0], vec![0.0, 1.0, 0.0], vec![0.0, 0.0, 1.0]];
        let r = Executor::new(prob, BFGS::new(MoreThuenteLineSearch::new()))
            .configure(|s| s.param(start.clone()).inv_hessian(eye).max_iters(100)).run()?;
        let p = r.state().get_best_param().unwrap();
        println!("{:<7} BFGS        r {:.6}  K {:.4}  SSE {:.6}", unit, p[0], p[1], r.state().get_best_cost());
    }
    Ok(())
}
