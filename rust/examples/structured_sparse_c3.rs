// C3: the Thomas algorithm (hand-written) for the morphogen gradient,
// checked by residual, mass balance, decay length and nalgebra's dense LU.
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
    let x4 = thomas(&[0.0, -1.0, -1.0, -1.0], &[2.0, 3.0, 3.0, 2.0],
                    &[-1.0, -1.0, -1.0, 0.0], &[1.0, 0.0, 0.0, 0.0]);
    println!("{:?}", x4.iter().map(|v| v * 21.0).collect::<Vec<_>>()); // B2

    let (n, k) = (100, 0.01);
    let (a, b, c, s) = morphogen(n, 1.0, k);
    let x = thomas(&a, &b, &c, &s);
    let mut res: f64 = 0.0;                              // residual max |A x - s|
    for i in 0..n {
        let mut ax = b[i] * x[i];
        if i > 0 { ax += a[i] * x[i - 1]; }
        if i < n - 1 { ax += c[i] * x[i + 1]; }
        res = res.max((ax - s[i]).abs());
    }
    println!("{:e}", res);
    println!("{}", k * x.iter().sum::<f64>());           // mass balance
    println!("{}", x[10] / x[0]);                        // one decay length
    println!("{}", (-1.0_f64).exp());
    let am = DMatrix::from_fn(n, n, |i, j| {
        if i == j { b[i] } else if i == j + 1 { a[i] } else if j == i + 1 { c[i] } else { 0.0 }
    });
    let xs = am.lu().solve(&DVector::from_vec(s.clone())).unwrap();
    let diff = x.iter().zip(xs.iter()).map(|(p, q)| (p - q).abs()).fold(0.0, f64::max);
    println!("{:e}", diff);
}
