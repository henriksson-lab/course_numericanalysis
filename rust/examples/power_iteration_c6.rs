// C6: inverse iteration with a shift; LU factorisation computed once.
use nalgebra::{DMatrix, DVector};

fn inv_it(a: &DMatrix<f64>, sigma: f64, x0: &[f64], tol: f64, maxit: usize) -> (f64, DVector<f64>, usize) {
    let n = a.nrows();
    let lu = (a - sigma * DMatrix::<f64>::identity(n, n)).lu();   // factor once
    let mut x = DVector::from_row_slice(x0);
    let (mut lambda, mut k) = (0.0, 0);
    while k < maxit {
        k += 1;
        let y = lu.solve(&x).expect("singular matrix");             // cheap solve each step
        x = &y / y.norm();
        let ax = a * &x;
        lambda = x.dot(&ax);
        if (&ax - lambda * &x).norm() < tol {
            break;
        }
    }
    (lambda, x, k)
}

fn main() {
    let p = DMatrix::from_row_slice(3, 3, &[
        0.80, 0.15, 0.05,
        0.05, 0.90, 0.05,
        0.05, 0.10, 0.85,
    ]);
    let id = DMatrix::<f64>::identity(3, 3);
    let pslow = &id + 0.01 * (&p - &id);
    let (_, x, k) = inv_it(&pslow.transpose(), 1.0 + 1e-6, &[1.0, 0.0, 0.0], 1e-12, 100);
    println!("{} {:?}", k, (&x / x.sum()).as_slice());
    for s in [0.82, 0.77] {
        let (lam, x, k) = inv_it(&p.transpose(), s, &[1.0, 0.0, 0.0], 1e-12, 100);
        let xr: Vec<f64> = x.iter().map(|v| (v * 1e4).round() / 1e4).collect();
        println!("{} {} {} {:?}", s, k, lam, xr);
    }
    println!("{}", 0.01f64.ln() / 0.8f64.ln());
}
