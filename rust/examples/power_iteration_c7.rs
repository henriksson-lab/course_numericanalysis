// C7: Rayleigh quotient iteration on the trait correlation matrix.
use nalgebra::{DMatrix, DVector};

fn rqi(a: &DMatrix<f64>, x0: &[f64], tol: f64, maxit: usize) {
    let n = a.nrows();
    let mut x = DVector::from_row_slice(x0);
    x /= x.norm();
    for k in 0..=maxit {
        let ax = a * &x;
        let sigma = x.dot(&ax);                       // Rayleigh quotient = new shift
        let r = (&ax - sigma * &x).norm();
        println!("{} {} {:.2e}", k, sigma, r);
        if r < tol {                                  // stop before A - sigma I is singular
            break;
        }
        let y = (a - sigma * DMatrix::<f64>::identity(n, n)).lu().solve(&x).unwrap();
        x = &y / y.norm();
    }
}

fn main() {
    let c = DMatrix::from_row_slice(3, 3, &[1.0, 0.4, 0.3, 0.4, 1.0, 0.2, 0.3, 0.2, 1.0]);
    rqi(&c, &[1.0, 1.0, 1.0], 1e-9, 20);
    rqi(&c, &[1.0, 0.0, 0.0], 1e-9, 20);
    rqi(&c, &[1.0, -1.0, 0.0], 1e-9, 20);
}
