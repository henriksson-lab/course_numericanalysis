// rust/examples/gradient_methods_c3.rs  —  cargo run --release --example gradient_methods_c3
// C3: steepest descent on the 4PL; the Hessian as a finite-difference Jacobian
// of the gradient; its eigenvalues (nalgebra) and condition number.
use nalgebra::{DMatrix, SymmetricEigen};
include!("gradient_methods_shared/defs.rs");

fn sd_armijo(mut x: Vec<f64>, tol: f64, maxit: usize) -> (Vec<f64>, usize) {
    let mut k = 0;
    while k < maxit {
        let g = grad_sse(&x);
        if norm(&g) < tol { break; }
        let (fx, g2) = (sse(&x), dot(&g, &g));
        let mut a = 1.0;
        while !(sse(&axpy(&x, -a, &g)) <= fx - 1e-4 * a * g2) { a /= 2.0; }
        x = axpy(&x, -a, &g);
        k += 1;
    }
    (x, k)
}
fn fd_jacobian(f: impl Fn(&[f64]) -> Vec<f64>, x: &[f64]) -> DMatrix<f64> {
    let f0 = f(x);
    let mut j = DMatrix::zeros(f0.len(), x.len());
    for c in 0..x.len() {
        let h = f64::EPSILON.sqrt() * x[c].abs().max(1.0);
        let mut xh = x.to_vec(); xh[c] += h;
        let fh = f(&xh);
        for r in 0..f0.len() { j[(r, c)] = (fh[r] - f0[r]) / h; }
    }
    j
}

fn main() {
    let (x, steps) = sd_armijo(P0.to_vec(), 1e-4, 10_000);
    println!("{}  {:.6?}  {:.4}", steps, x, sse(&x));
    let h = fd_jacobian(grad_sse, &x);
    println!("{:.1}", h);
    println!("{:.2e}", (&h - h.transpose()).abs().max());
    let mut ev: Vec<f64> = SymmetricEigen::new((&h + h.transpose()) / 2.0).eigenvalues.iter().cloned().collect();
    ev.sort_by(|a, b| b.partial_cmp(a).unwrap());
    let kap = ev[0] / ev[3];
    println!("{:.6?}  kappa {:.2}", ev, kap);
    let q = (kap - 1.0) / (kap + 1.0);
    println!("{:.7}  {:.0}", q, 1e-6f64.ln() / q.ln());
}
