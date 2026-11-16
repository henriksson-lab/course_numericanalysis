// C1: check the hand-derived Jacobian against a finite-difference Jacobian.
use nalgebra::{DMatrix, DVector};

const ALPHA: f64 = 3.0;
const N: f64 = 2.0;

fn f(x: &DVector<f64>) -> DVector<f64> {
    let (u, v) = (x[0], x[1]);
    DVector::from_vec(vec![ALPHA / (1.0 + v.powf(N)) - u,
                           ALPHA / (1.0 + u.powf(N)) - v])
}

fn jac(x: &DVector<f64>) -> DMatrix<f64> {
    let (u, v) = (x[0], x[1]);
    DMatrix::from_row_slice(2, 2, &[
        -1.0, -ALPHA * N * v.powf(N - 1.0) / (1.0 + v.powf(N)).powi(2),
        -ALPHA * N * u.powf(N - 1.0) / (1.0 + u.powf(N)).powi(2), -1.0])
}

// works for any F: R^n -> R^m (here passed as a function)
fn fd_jacobian(f: impl Fn(&DVector<f64>) -> DVector<f64>, x: &DVector<f64>, h: f64) -> DMatrix<f64> {
    let fx = f(x);
    let mut jm = DMatrix::zeros(fx.len(), x.len());
    for j in 0..x.len() {
        let hj = h * x[j].abs().max(1.0);       // step relative to the size of x[j]
        let mut xh = x.clone();
        xh[j] += hj;                             // nudge only x[j]
        jm.set_column(j, &((f(&xh) - &fx) / hj)); // column j
    }
    jm
}

fn main() {
    let h = f64::EPSILON.sqrt();
    println!("{}", jac(&DVector::from_vec(vec![1.0, 2.0]))); // B2
    for p in [[3.0, 0.5], [1.0, 2.0], [1.2, 1.2], [0.1, 4.0]] {
        let x = DVector::from_vec(p.to_vec());
        let diff = (fd_jacobian(f, &x, h) - jac(&x)).abs().max();
        println!("{:?} : {:e}", p, diff);
    }
}
