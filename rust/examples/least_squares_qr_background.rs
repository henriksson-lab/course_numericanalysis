// Background: a straight-line fit by the normal equations and by QR (nalgebra).
use nalgebra::{DMatrix, DVector};

fn main() {
    let log_n: Vec<f64> = (1..=7).flat_map(|k| [k as f64; 3]).collect();
    let cq = DVector::from_vec(vec![
        34.52, 34.69, 34.48, 31.62, 31.37, 31.14, 28.05, 28.10, 28.07,
        24.54, 24.90, 24.68, 21.13, 20.81, 21.47, 17.89, 17.90, 18.09,
        14.71, 14.67, 14.73]);
    // design matrix, 21 x 2: a column of ones and log10 N
    let a = DMatrix::from_fn(log_n.len(), 2, |i, j| if j == 0 { 1.0 } else { log_n[i] });

    // normal equations (A^T A) x = A^T b, solved by Cholesky
    let x_ne = (a.transpose() * &a).cholesky().unwrap().solve(&(a.transpose() * &cq));
    println!("{:.6} {:.6}", x_ne[0], x_ne[1]);

    // Householder QR: R x = Q^T b; Q is never formed
    let qr = a.clone().qr();
    let mut qtb = cq.clone();
    qr.q_tr_mul(&mut qtb);                                    // qtb <- Q^T b
    let x_qr = qr.r().solve_upper_triangular(&qtb.rows(0, 2).into_owned()).unwrap();
    println!("{:.6} {:.6}", x_qr[0], x_qr[1]);

    // 2-norm condition number = largest / smallest singular value
    let sv = a.singular_values();
    println!("{}", sv.max() / sv.min());

    // BCA standards: BSA in ug/mL (duplicates) and absorbance at 562 nm
    let conc: Vec<f64> = [0., 25., 125., 250., 500., 750., 1000., 1500., 2000.]
        .iter().flat_map(|&c| [c; 2]).collect();
    let a562 = DVector::from_vec(vec![0.087, 0.103, 0.150, 0.109, 0.226, 0.229, 0.358,
        0.344, 0.601, 0.569, 0.780, 0.789, 0.951, 0.942, 1.300, 1.238, 1.545, 1.532]);
    println!("{} {}", conc.len(), a562.len());
}
