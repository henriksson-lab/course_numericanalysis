// C1: power iteration for the fish Leslie matrix and the 2 x 2 plant.
use nalgebra::{DMatrix, DVector};

// returns (lambda, x, iterations, residual)
fn power_it(a: &DMatrix<f64>, x0: &[f64], tol: f64, maxit: usize) -> (f64, DVector<f64>, usize, f64) {
    let mut x = DVector::from_row_slice(x0);
    x /= x.sum();
    let (mut lambda, mut r, mut k) = (0.0, f64::INFINITY, 0);
    while k < maxit {
        k += 1;
        let y = a * &x;
        lambda = y.sum();                       // growth factor of the total
        r = (&y - lambda * &x).norm();          // residual ||Ax - lambda x||
        x = y / lambda;                         // normalise: proportions
        if r < tol {
            break;
        }
    }
    (lambda, x, k, r)
}

fn main() {
    let l = DMatrix::from_row_slice(4, 4, &[
        0.0, 1.5, 4.0, 5.0,
        0.3, 0.0, 0.0, 0.0,
        0.0, 0.6, 0.0, 0.0,
        0.0, 0.0, 0.5, 0.0,
    ]);
    let (lam, x, k, r) = power_it(&l, &[0.0, 1000.0, 0.0, 0.0], 1e-10, 1000);
    println!("{} {:.8?} {} {:e}", lam, x.as_slice(), k, r);
    // compare with the largest eigenvalue from nalgebra (eigenvalues only)
    let lam1 = l.complex_eigenvalues().iter().map(|z| z.re).fold(f64::MIN, f64::max);
    println!("{:e}", (lam - lam1).abs());

    let plant = DMatrix::from_row_slice(2, 2, &[0.6, 0.9, 0.8, 0.0]);
    let (lam, x, k, _) = power_it(&plant, &[0.0, 1.0], 1e-10, 1000);
    println!("{} {:?} {}", lam, x.as_slice(), k);
}
