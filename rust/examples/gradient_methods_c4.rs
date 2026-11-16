// rust/examples/gradient_methods_c4.rs  —  cargo run --example gradient_methods_c4
// C4: Newton for minimisation with backtracking; a start with an indefinite
// Hessian; modified Newton (H + mu I). Linear algebra with nalgebra.
use nalgebra::{DMatrix, DVector, SymmetricEigen};
include!("gradient_methods_shared/defs.rs");

fn hess(p: &[f64]) -> DMatrix<f64> {   // central-difference Jacobian of the gradient
    let mut h = DMatrix::zeros(4, 4);
    for j in 0..4 {
        let s = 1e-5 * p[j].abs().max(1.0);
        let (mut a, mut b) = (p.to_vec(), p.to_vec());
        a[j] += s; b[j] -= s;
        let (ga, gb) = (grad_sse(&a), grad_sse(&b));
        for i in 0..4 { h[(i, j)] = (ga[i] - gb[i]) / (2.0 * s); }
    }
    h
}
fn eig_sorted(h: &DMatrix<f64>) -> Vec<f64> {
    let mut ev: Vec<f64> = SymmetricEigen::new((h + h.transpose()) / 2.0).eigenvalues.iter().cloned().collect();
    ev.sort_by(|a, b| b.partial_cmp(a).unwrap());
    ev
}
fn newton_dir(h: DMatrix<f64>, g: &[f64]) -> Vec<f64> {
    let s = h.lu().solve(&(-DVector::from_column_slice(g))).expect("singular Hessian");
    s.iter().cloned().collect()
}
fn backtrack(x: &[f64], g: &[f64], s: &[f64]) -> f64 {
    let (fx, gs) = (sse(x), dot(g, s));
    let mut a = 1.0;
    while !(sse(&axpy(x, a, s)) <= fx + 1e-4 * a * gs) { a /= 2.0; }
    a
}

fn main() {
    // (a)
    let mut x = P0.to_vec();
    for k in 0..6 {
        let g = grad_sse(&x);
        println!("{}  SSE {:.8}  |g| {:.2e}", k, sse(&x), norm(&g));
        if norm(&g) < 1e-4 { break; }
        let s = newton_dir(hess(&x), &g);
        x = axpy(&x, backtrack(&x, &g, &s), &s);
    }
    // (b)
    let mut x = vec![0.0, 100.0, 10.0, 1.0];
    println!("{:.5?}", eig_sorted(&hess(&x)));
    for _ in 0..20 { let s = newton_dir(hess(&x), &grad_sse(&x)); x = axpy(&x, 1.0, &s); }
    println!("{:.2?}  {:.2}", x, sse(&x));
    // (c)
    let mut x = vec![0.0, 100.0, 10.0, 1.0];
    let mut steps = 0;
    while steps < 100 {
        let g = grad_sse(&x);
        if norm(&g) < 1e-4 { break; }
        let mut h = hess(&x);
        h = (&h + h.transpose()) / 2.0;
        let lmin = *eig_sorted(&h).last().unwrap();
        if lmin <= 0.0 {
            let mu = 1e-3 * h.diagonal().abs().max() - lmin;
            h += DMatrix::identity(4, 4) * mu;
        }
        let s = newton_dir(h, &g);
        x = axpy(&x, backtrack(&x, &g, &s), &s);
        steps += 1;
    }
    println!("{}  {:.6?}  {:.4}", steps, x, sse(&x));
}
