// C3: power iteration fails when |lambda2| = |lambda1|; and why we normalise.
use nalgebra::{DMatrix, DVector};

fn power_it(a: &DMatrix<f64>, x0: &[f64], tol: f64, maxit: usize) -> (f64, usize, f64) {
    let mut x = DVector::from_row_slice(x0);
    x /= x.sum();
    let (mut lambda, mut r, mut k) = (0.0, f64::INFINITY, 0);
    while k < maxit {
        k += 1;
        let y = a * &x;
        lambda = y.sum();
        r = (&y - lambda * &x).norm();
        x = y / lambda;
        if r < tol {
            break;
        }
    }
    (lambda, k, r)
}

fn main() {
    let ls = DMatrix::from_row_slice(4, 4, &[
        0.0, 0.0, 0.0, 12.0,                         // spawn once, at age 3
        0.3, 0.0, 0.0, 0.0,
        0.0, 0.6, 0.0, 0.0,
        0.0, 0.0, 0.5, 0.0,
    ]);
    let m: Vec<f64> = ls.complex_eigenvalues().iter().map(|z| z.norm()).collect();
    println!("{:?}", m);
    let (_, k, r) = power_it(&ls, &[0.0, 1.0, 0.0, 0.0], 1e-10, 200);
    println!("{} {}", k, r);
    let biennial = DMatrix::from_row_slice(2, 2, &[0.0, 4.0, 0.25, 0.0]);
    println!("{:?}", power_it(&biennial, &[1.0, 0.0], 1e-10, 200));
    let early = DMatrix::from_row_slice(2, 2, &[0.4, 4.0, 0.25, 0.0]);
    println!("{:?}", early.complex_eigenvalues().iter().map(|z| z.re).collect::<Vec<f64>>());
    println!("{:?}", power_it(&early, &[1.0, 0.0], 1e-10, 1000));

    let l = DMatrix::from_row_slice(4, 4, &[
        0.0, 1.5, 4.0, 5.0,
        0.3, 0.0, 0.0, 0.0,
        0.0, 0.6, 0.0, 0.0,
        0.0, 0.0, 0.5, 0.0,
    ]);
    let mut x = DVector::from_row_slice(&[0.0, 1.0, 0.0, 0.0]);
    for _ in 0..5000 {
        x = &l * &x;                                 // no normalisation
    }
    println!("{:?}", x.as_slice());
}
