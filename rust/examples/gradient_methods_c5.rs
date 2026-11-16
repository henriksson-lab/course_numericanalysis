// rust/examples/gradient_methods_c5.rs  —  cargo run --example gradient_methods_c5
// C5: BFGS written out, with call counters, against argmin's BFGS.
use argmin::core::{Executor, State};
use argmin::solver::linesearch::MoreThuenteLineSearch;
use argmin::solver::quasinewton::BFGS;
use std::cell::Cell;
use std::sync::atomic::{AtomicU64, Ordering};
include!("gradient_methods_shared/defs.rs");

// argmin owns the problem, so count its calls with global counters
static NF: AtomicU64 = AtomicU64::new(0);
static NG: AtomicU64 = AtomicU64::new(0);
struct Counted;
impl CostFunction for Counted {
    type Param = Vec<f64>;
    type Output = f64;
    fn cost(&self, p: &Vec<f64>) -> Result<f64, Error> { NF.fetch_add(1, Ordering::Relaxed); Ok(sse(p)) }
}
impl Gradient for Counted {
    type Param = Vec<f64>;
    type Gradient = Vec<f64>;
    fn gradient(&self, p: &Vec<f64>) -> Result<Vec<f64>, Error> { NG.fetch_add(1, Ordering::Relaxed); Ok(grad_sse(p)) }
}

fn bfgs(f: &dyn Fn(&[f64]) -> f64, g: &dyn Fn(&[f64]) -> Vec<f64>,
        mut x: Vec<f64>, tol: f64, maxit: usize) -> (Vec<f64>, usize) {
    let n = x.len();
    let mut b = vec![vec![0.0; n]; n];          // B ~ inverse Hessian, start = I
    for i in 0..n { b[i][i] = 1.0; }
    let mut gk = g(&x);
    let mut k = 0;
    while k < maxit && norm(&gk) >= tol {
        let s: Vec<f64> = (0..n).map(|i| -dot(&b[i], &gk)).collect();
        let (fx, gs) = (f(&x), dot(&gk, &s));
        let mut a = 1.0;
        while !(f(&axpy(&x, a, &s)) <= fx + 1e-4 * a * gs) { a /= 2.0; }
        let s: Vec<f64> = s.iter().map(|v| a * v).collect();
        x = axpy(&x, 1.0, &s);
        let gn = g(&x);
        let yv: Vec<f64> = gn.iter().zip(&gk).map(|(a, b)| a - b).collect();
        gk = gn;
        let sy = dot(&yv, &s);
        if sy > 0.0 {   // B <- V^T B V + s s^T / sy,  V = I - y s^T / sy
            let v: Vec<Vec<f64>> = (0..n).map(|i| (0..n).map(|j|
                (i == j) as i32 as f64 - yv[i] * s[j] / sy).collect()).collect();
            let bv: Vec<Vec<f64>> = (0..n).map(|i| (0..n).map(|j|
                (0..n).map(|l| b[i][l] * v[l][j]).sum()).collect()).collect();
            b = (0..n).map(|i| (0..n).map(|j|
                (0..n).map(|l| v[l][i] * bv[l][j]).sum::<f64>() + s[i] * s[j] / sy).collect()).collect();
        }
        k += 1;
    }
    (x, k)
}

fn main() -> Result<(), Error> {
    let (nf, ng) = (Cell::new(0), Cell::new(0));
    let fc = |p: &[f64]| { nf.set(nf.get() + 1); sse(p) };
    let gc = |p: &[f64]| { ng.set(ng.get() + 1); grad_sse(p) };
    let (x, steps) = bfgs(&fc, &gc, P0.to_vec(), 1e-4, 200);
    println!("own BFGS: steps {}  sse calls {}  gradient calls {}  SSE {:.6}", steps, nf.get(), ng.get(), sse(&x));
    let eye: Vec<Vec<f64>> = (0..4).map(|i| (0..4).map(|j| (i == j) as i32 as f64).collect()).collect();
    let res = Executor::new(Counted, BFGS::new(MoreThuenteLineSearch::new()))
        .configure(|s| s.param(P0.to_vec()).inv_hessian(eye).max_iters(200))
        .run()?;
    let s = res.state();
    println!("argmin:   iterations {}  sse calls {}  gradient calls {}  SSE {:.6}  |g| {:.1e}",
             s.get_iter(), NF.load(Ordering::Relaxed), NG.load(Ordering::Relaxed), s.get_best_cost(),
             norm(&grad_sse(s.get_best_param().unwrap())));
    Ok(())
}
