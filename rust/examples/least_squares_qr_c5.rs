// C5: error against kappa, normal equations vs Householder QR.
use nalgebra::{DMatrix, DVector};

// scientific notation as R and Python print it: 2.2e+01 rather than Rust's 2.2e1
fn sci(x: f64, digits: usize) -> String {
    let s = format!("{:.*e}", digits, x);
    let (m, e) = s.split_once('e').unwrap();
    let e: i32 = e.parse().unwrap();
    format!("{}e{}{:02}", m, if e < 0 { '-' } else { '+' }, e.abs())
}

fn main() {
    let t: Vec<f64> = (0..50).map(|i| i as f64 / 49.0).collect();
    for d in 1..=12usize {
        let a = DMatrix::from_fn(50, d + 1, |i, j| t[i].powi(j as i32));
        let xt = DVector::from_element(d + 1, 1.0);
        let b = &a * &xt;                                      // data with known coefficients
        let relerr = |x: &DVector<f64>| (x - &xt).norm() / xt.norm();
        let ne = (a.transpose() * &a).lu().solve(&(a.transpose() * &b)).unwrap();
        let qr = a.clone().qr();
        let mut qtb = b.clone();
        qr.q_tr_mul(&mut qtb);
        let xq = qr.r().solve_upper_triangular(&qtb.rows(0, d + 1).into_owned()).unwrap();
        let sv = a.singular_values();
        println!("{:2} {:>8} {:>8} {:>8}", d, sci(sv.max() / sv.min(), 1),
                 sci(relerr(&ne), 1), sci(relerr(&xq), 1));
    }
}
