// C4: timing Thomas against a dense LU solve (nalgebra). Run with --release.
use nalgebra::{DMatrix, DVector};
use std::time::Instant;

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

/// Morphogen in a row of n cells: exchange D, decay k, production 1 in cell 1.
fn morphogen(n: usize, dd: f64, k: f64) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let mut a = vec![-dd; n];
    a[0] = 0.0;
    let mut c = vec![-dd; n];
    c[n - 1] = 0.0;
    let mut b = vec![2.0 * dd + k; n];
    b[0] = dd + k;
    b[n - 1] = dd + k;
    let mut s = vec![0.0; n];
    s[0] = 1.0;
    (a, b, c, s)
}

fn main() {
    for n in [1000, 2000, 4000] {
        let (a, b, c, s) = morphogen(n, 1.0, 0.01);
        let am = DMatrix::from_fn(n, n, |i, j| {
            if i == j { b[i] } else if i == j + 1 { a[i] } else if j == i + 1 { c[i] } else { 0.0 }
        });
        let t0 = Instant::now();
        let _ = am.lu().solve(&DVector::from_vec(s.clone()));
        let t_dense = t0.elapsed().as_secs_f64();
        let t1 = Instant::now();
        let _ = thomas(&a, &b, &c, &s);
        println!("{} {:.3} {:.6}", n, t_dense, t1.elapsed().as_secs_f64());
    }
    for n in [100_000, 1_000_000] {
        let (a, b, c, s) = morphogen(n, 1.0, 0.01);
        let t1 = Instant::now();
        let _ = thomas(&a, &b, &c, &s);
        println!("{} {:.6}", n, t1.elapsed().as_secs_f64());
    }
}
