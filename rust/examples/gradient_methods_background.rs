// Background: the 4PL dose-response fit. Data, SSE, analytic gradient, a
// central-difference check, and argmin's BFGS (the counterpart of
// optim(method = "BFGS") / scipy's minimize(method="BFGS")).
use argmin::core::{CostFunction, Error, Executor, Gradient, State};
use argmin::solver::linesearch::MoreThuenteLineSearch;
use argmin::solver::quasinewton::BFGS;

const Y: [f64; 18] = [12.9, 9.8, 2.9, 0.1, 4.2, 9.7, 17.7, 22.3, 49.7, 42.1,
                      67.3, 72.6, 97.4, 91.8, 94.6, 93.3, 92.9, 95.1];

fn doses() -> Vec<f64> {          // 9 doses 0.01 ... 100 µM, 2 wells each
    (0..18).map(|i| 10f64.powf(-2.0 + 0.5 * (i / 2) as f64)).collect()
}
fn f4(p: &[f64], c: f64) -> f64 {  // p = (bottom, top, EC50, hill)
    p[0] + (p[1] - p[0]) / (1.0 + (p[2] / c).powf(p[3]))
}
fn sse(p: &[f64]) -> f64 {
    doses().iter().zip(Y.iter()).map(|(&c, &y)| (y - f4(p, c)).powi(2)).sum()
}
fn grad_sse(p: &[f64]) -> Vec<f64> {
    let mut g = vec![0.0; 4];
    for (&c, &y) in doses().iter().zip(Y.iter()) {
        let u = (p[2] / c).powf(p[3]);
        let d = 1.0 + u;
        let r = y - f4(p, c);
        let j = [u / d, 1.0 / d,
                 -(p[1] - p[0]) * p[3] * u / (p[2] * d * d),
                 -(p[1] - p[0]) * u * (p[2] / c).ln() / (d * d)];
        for k in 0..4 { g[k] -= 2.0 * j[k] * r; }
    }
    g
}

struct Fit;   // argmin wants the problem as a type with cost and gradient
impl CostFunction for Fit {
    type Param = Vec<f64>;
    type Output = f64;
    fn cost(&self, p: &Vec<f64>) -> Result<f64, Error> { Ok(sse(p)) }
}
impl Gradient for Fit {
    type Param = Vec<f64>;
    type Gradient = Vec<f64>;
    fn gradient(&self, p: &Vec<f64>) -> Result<Vec<f64>, Error> { Ok(grad_sse(p)) }
}

fn main() -> Result<(), Error> {
    let p0 = vec![0.0, 100.0, 1.0, 1.0];
    println!("{:.4}", sse(&p0));
    println!("{:?}", grad_sse(&p0));
    // central differences, step 1e-6 * max(|p_j|, 1)
    let fd: Vec<f64> = (0..4).map(|j| {
        let h = 1e-6 * p0[j].abs().max(1.0);
        let (mut a, mut b) = (p0.clone(), p0.clone());
        a[j] += h; b[j] -= h;
        (sse(&a) - sse(&b)) / (2.0 * h)
    }).collect();
    println!("{:?}", fd);
    let identity: Vec<Vec<f64>> = (0..4).map(|i| (0..4).map(|j| if i == j { 1.0 } else { 0.0 }).collect()).collect();
    let solver = BFGS::new(MoreThuenteLineSearch::new());
    let res = Executor::new(Fit, solver)
        .configure(|s| s.param(p0.clone()).inv_hessian(identity).max_iters(200))
        .run()?;
    let s = res.state();
    println!("{:?}  SSE {:.4}  iterations {}", s.get_best_param().unwrap(), s.get_best_cost(), s.get_iter());
    Ok(())
}
