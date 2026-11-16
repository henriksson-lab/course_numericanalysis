// C5: sum ratio against Rayleigh quotient, symmetric and nonsymmetric.
use nalgebra::{DMatrix, DVector};

// first iteration at which each estimate is within 1e-10 of lam1
fn compare(a: &DMatrix<f64>, x0: &[f64], lam1: f64, kmax: usize) -> (Option<usize>, Option<usize>) {
    let mut x = DVector::from_row_slice(x0);
    let (mut first_sum, mut first_rq) = (None, None);
    for k in 1..=kmax {
        let y = a * &x;
        let e_sum = (y.sum() / x.sum() - lam1).abs();
        let e_rq = (x.dot(&y) / x.dot(&x) - lam1).abs();
        if e_sum < 1e-10 && first_sum.is_none() { first_sum = Some(k); }
        if e_rq < 1e-10 && first_rq.is_none() { first_rq = Some(k); }
        x = &y / y.norm();
    }
    (first_sum, first_rq)
}

fn main() {
    let c = DMatrix::from_row_slice(3, 3, &[1.0, 0.4, 0.3, 0.4, 1.0, 0.2, 0.3, 0.2, 1.0]);
    let lam_c = c.clone().symmetric_eigen().eigenvalues.max();
    println!("{:?}", compare(&c, &[1.0, 0.0, 0.0], lam_c, 80));
    let l = DMatrix::from_row_slice(4, 4, &[
        0.0, 1.5, 4.0, 5.0,
        0.3, 0.0, 0.0, 0.0,
        0.0, 0.6, 0.0, 0.0,
        0.0, 0.0, 0.5, 0.0,
    ]);
    let lam_l = l.complex_eigenvalues().iter().map(|z| z.norm()).fold(0.0, f64::max);
    println!("{:?}", compare(&l, &[0.0, 1.0, 0.0, 0.0], lam_l, 80));
}
