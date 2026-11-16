// C2: reparameterisation with exp and softmax, minimised with argmin's Nelder-Mead
// (R and SciPy use BFGS, so the last digits of the parameters differ).
use argmin::core::{CostFunction, Error, Executor, State};
use argmin::solver::neldermead::NelderMead;

struct Decay; // (a) SSE of A e^{-kt} + exp(theta)
impl CostFunction for Decay {
    type Param = Vec<f64>;
    type Output = f64;
    fn cost(&self, q: &Vec<f64>) -> Result<f64, Error> {
        let y = [1.020, 0.577, 0.378, 0.243, 0.125, 0.062, 0.020, 0.010, -0.012];
        Ok(y.iter().enumerate()
            .map(|(i, yi)| (yi - q[0] * (-q[1] * i as f64).exp() - q[2].exp()).powi(2))
            .sum())
    }
}

fn softmax(th: &[f64]) -> [f64; 3] {
    let e = [1.0, th[0].exp(), th[1].exp()]; // first theta fixed at 0
    let s: f64 = e.iter().sum();
    [e[0] / s, e[1] / s, e[2] / s]
}

struct Diet; // (b) isotope misfit of the diet fractions softmax(0, th2, th3)
impl CostFunction for Diet {
    type Param = Vec<f64>;
    type Output = f64;
    fn cost(&self, th: &Vec<f64>) -> Result<f64, Error> {
        let a = [[-26.0, -10.0, -18.0], [2.0, 6.0, 14.0]]; // d13C, d15N
        let obs = [-24.0, 9.0];
        let p = softmax(th);
        Ok((0..2).map(|i| ((0..3).map(|j| a[i][j] * p[j]).sum::<f64>() - obs[i]).powi(2)).sum())
    }
}

fn simplex(start: &[f64]) -> Vec<Vec<f64>> {
    let mut s = vec![start.to_vec()];
    for i in 0..start.len() {
        let mut v = start.to_vec();
        v[i] += 0.5;
        s.push(v);
    }
    s
}

fn main() -> Result<(), Error> {
    let solver = NelderMead::new(simplex(&[1.0, 0.3, 0.1f64.ln()])).with_sd_tolerance(1e-14)?;
    let res = Executor::new(Decay, solver).configure(|s| s.max_iters(2000)).run()?;
    let q = res.state().get_best_param().unwrap().clone();
    println!("theta: {:.3} {:.3} {:.3}  c = {:.3e}  iterations {}",
             q[0], q[1], q[2], q[2].exp(), res.state().get_iter());

    let solver = NelderMead::new(simplex(&[0.0, 0.0])).with_sd_tolerance(1e-14)?;
    let res = Executor::new(Diet, solver).configure(|s| s.max_iters(2000)).run()?;
    let th = res.state().get_best_param().unwrap().clone();
    let p = softmax(&th);
    println!("p = {:.4} {:.4} {:.4}  loss {:.4}  theta {:.2} {:.2}",
             p[0], p[1], p[2], res.state().get_best_cost(), th[0], th[1]);
    let pe = [27.0 / 52.0, 0.0, 25.0 / 52.0]; // exact answer (maize = 0)
    let r = [-26.0 * pe[0] - 18.0 * pe[2] + 24.0, 2.0 * pe[0] + 14.0 * pe[2] - 9.0];
    println!("exact loss {:.4}", r[0] * r[0] + r[1] * r[1]);
    Ok(())
}
