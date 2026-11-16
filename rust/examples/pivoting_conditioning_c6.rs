// Residual vs error for an ill-conditioned 6-dye panel (peaks 6 nm apart).
use nalgebra::{DMatrix, DVector};

fn panel(n: usize, spacing: f64, width: f64) -> DMatrix<f64> {
    let peak = |i: usize| 500.0 + spacing * i as f64;
    DMatrix::from_fn(n, n, |i, j| (-(peak(i) - peak(j)).powi(2) / (2.0 * width * width)).exp())
}

fn main() {
    let xt = DVector::from_vec(vec![5.0, 1.0, 3.0, 2.0, 4.0, 1.0]);
    let a = panel(6, 6.0, 30.0);
    let b = &a * &xt;
    let xh = a.clone().lu().solve(&b).unwrap();
    let sv = a.singular_values();
    let k = sv.max() / sv.min();
    let r = (&b - &a * &xh).norm() / b.norm(); // relative residual
    let e = (&xh - &xt).norm() / xt.norm();     // relative error
    println!("kappa = {:.2e}  residual = {:.2e}  error = {:.2e}  bound = {:.2e}", k, r, e, k * r);
}

// Output:
// kappa = 1.14e8  residual = 8.77e-17  error = 1.18e-9  bound = 1.00e-8
