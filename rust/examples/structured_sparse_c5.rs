// C5: Thomas has no pivoting. A tiny first pivot gives a wrong answer;
// the residual exposes it. nalgebra's LU (partial pivoting) gets it right.
use nalgebra::{DMatrix, DVector};

/// Thomas algorithm: a = below the diagonal (a[0] unused), b = diagonal,
/// c = above the diagonal (c[n-1] unused), d = right-hand side.
fn thomas(a: &[f64], b: &[f64], c: &[f64], d: &[f64]) -> Vec<f64> {
    let n = d.len();
    let (mut b, mut d) = (b.to_vec(), d.to_vec()); // copies: inputs stay unchanged
    for i in 1..n {
        let w = a[i] / b[i - 1];
        b[i] -= w * c[i - 1];
        d[i] -= w * d[i - 1];
    }
    let mut x = vec![0.0; n];
    x[n - 1] = d[n - 1] / b[n - 1];
    for i in (0..n - 1).rev() {
        x[i] = (d[i] - c[i] * x[i + 1]) / b[i];
    }
    x
}

fn diag_dominant(a: &[f64], b: &[f64], c: &[f64]) -> bool {
    (0..b.len()).all(|i| b[i].abs() >= a[i].abs() + c[i].abs())
}

fn main() {
    let a = [0.0, 1.0, 1.0, 1.0];
    let b = [1e-17, 2.0, 2.0, 2.0];
    let c = [1.0, 1.0, 1.0, 0.0];
    let d = [1.0, 1.0, 1.0, 1.0];
    println!("{}", diag_dominant(&a, &b, &c));
    let am = DMatrix::from_fn(4, 4, |i, j| {
        if i == j { b[i] } else if i == j + 1 { a[i] } else if j == i + 1 { c[i] } else { 0.0 }
    });
    let x = thomas(&a, &b, &c, &d);
    println!("{:?}", x);
    let dv = DVector::from_row_slice(&d);
    println!("{:?}", am.clone().lu().solve(&dv).unwrap().as_slice()); // with pivoting
    let res = (&am * DVector::from_vec(x) - dv).abs().max();
    println!("{}", res);                                  // residual exposes it

    let (n, k) = (100, 0.01);                             // the morphogen matrix
    let mut a = vec![-1.0; n]; a[0] = 0.0;
    let mut c = vec![-1.0; n]; c[n - 1] = 0.0;
    let mut b = vec![2.0 + k; n]; b[0] = 1.0 + k; b[n - 1] = 1.0 + k;
    println!("{}", diag_dominant(&a, &b, &c));
}
