// rust/examples/gradient_methods_c1.rs  —  cargo run --example gradient_methods_c1
// C1: check the gradient; a missing factor; a sign error fed to argmin's BFGS.
use argmin::core::{Executor, State};
use argmin::solver::linesearch::MoreThuenteLineSearch;
use argmin::solver::quasinewton::BFGS;
include!("gradient_methods_shared/defs.rs");

fn fd_grad(f: impl Fn(&[f64]) -> f64, x: &[f64], central: bool) -> Vec<f64> {
    (0..x.len()).map(|j| {
        let (mut a, mut b) = (x.to_vec(), x.to_vec());
        if central {
            let h = 1e-5 * x[j].abs().max(1.0);
            a[j] += h; b[j] -= h;
            (f(&a) - f(&b)) / (2.0 * h)
        } else {
            let h = f64::EPSILON.sqrt() * x[j].abs().max(1.0);
            a[j] += h;
            (f(&a) - f(x)) / h
        }
    }).collect()
}
fn rel(a: &[f64], b: &[f64]) -> String {
    a.iter().zip(b).map(|(x, y)| format!("{:.1e}", ((x - y) / y).abs())).collect::<Vec<_>>().join("  ")
}
// the model Jacobian with two possible bugs
fn grad_buggy(p: &[f64], forget_hill: bool, flip_sign: bool) -> Vec<f64> {
    let mut g = vec![0.0; 4];
    for (&c, &y) in doses().iter().zip(Y.iter()) {
        let u = (p[2] / c).powf(p[3]);
        let d = 1.0 + u;
        let r = y - f4(p, c);
        let h = if forget_hill { 1.0 } else { p[3] };
        let l = if flip_sign { (c / p[2]).ln() } else { (p[2] / c).ln() };
        let j = [u / d, 1.0 / d, -(p[1] - p[0]) * h * u / (p[2] * d * d), -(p[1] - p[0]) * u * l / (d * d)];
        for k in 0..4 { g[k] -= 2.0 * j[k] * r; }
    }
    g
}
struct BadFit;
impl CostFunction for BadFit {
    type Param = Vec<f64>;
    type Output = f64;
    fn cost(&self, p: &Vec<f64>) -> Result<f64, Error> { Ok(sse(p)) }
}
impl Gradient for BadFit {
    type Param = Vec<f64>;
    type Gradient = Vec<f64>;
    fn gradient(&self, p: &Vec<f64>) -> Result<Vec<f64>, Error> { Ok(grad_buggy(p, false, true)) }
}

fn main() -> Result<(), Error> {
    let p0 = P0.to_vec();
    let g = grad_sse(&p0);
    // (a)
    println!("analytic {:.6?}", g);
    println!("forward  {:.6?}", fd_grad(sse, &p0, false));
    println!("{}", rel(&fd_grad(sse, &p0, false), &g));
    println!("{}", rel(&fd_grad(sse, &p0, true), &g));
    // (b) missing hill factor: invisible at hill = 1
    let p1 = [5.0, 95.0, 2.0, 1.5];
    println!("{}", rel(&grad_buggy(&p0, true, false), &fd_grad(sse, &p0, true)));
    println!("{}", rel(&grad_buggy(&p1, true, false), &fd_grad(sse, &p1, true)));
    // (c) BFGS with a sign error in dm/d hill
    let eye: Vec<Vec<f64>> = (0..4).map(|i| (0..4).map(|j| (i == j) as i32 as f64).collect()).collect();
    let res = Executor::new(BadFit, BFGS::new(MoreThuenteLineSearch::new()))
        .configure(|s| s.param(p0.clone()).inv_hessian(eye).max_iters(200))
        .run();
    match res {
        Ok(r) => { let s = r.state();
            println!("{:.4?}  SSE {:.4}  {}", s.get_best_param().unwrap(), s.get_best_cost(),
                     s.get_termination_reason().unwrap().text()); }
        Err(e) => println!("argmin stopped with an error: {}", e),
    }
    Ok(())
}
