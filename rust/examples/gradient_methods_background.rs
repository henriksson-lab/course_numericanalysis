// rust/examples/gradient_methods_background.rs  —  cargo run --example gradient_methods_background
// The 4PL fit: SSE, gradient, a central-difference check, and argmin's BFGS
// (the counterpart of optim(method = "BFGS") and minimize(method="BFGS")).
use argmin::core::{Executor, State};
use argmin::solver::linesearch::MoreThuenteLineSearch;
use argmin::solver::quasinewton::BFGS;
include!("gradient_methods_shared/defs.rs");

fn main() -> Result<(), Error> {
    let p0 = P0.to_vec();
    println!("{:.4}", sse(&p0));
    println!("{:.6?}", grad_sse(&p0));
    let fd: Vec<f64> = (0..4).map(|j| {           // central differences
        let h = 1e-6 * p0[j].abs().max(1.0);
        let (mut a, mut b) = (p0.clone(), p0.clone());
        a[j] += h; b[j] -= h;
        (sse(&a) - sse(&b)) / (2.0 * h)
    }).collect();
    println!("{:.6?}", fd);
    let eye: Vec<Vec<f64>> = (0..4).map(|i| (0..4).map(|j| (i == j) as i32 as f64).collect()).collect();
    let res = Executor::new(Fit, BFGS::new(MoreThuenteLineSearch::new()))
        .configure(|s| s.param(p0.clone()).inv_hessian(eye).max_iters(200))
        .run()?;
    let s = res.state();
    println!("{:.6?}  SSE {:.6}  iterations {}", s.get_best_param().unwrap(), s.get_best_cost(), s.get_iter());
    Ok(())
}
