// C2: Gauss-Newton with a forward-difference Jacobian and a QR least-squares solve
use nalgebra::{DMatrix, DVector};

const D: f64 = 100.0;
const TT: [f64; 12] = [0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 12.0];
const Y: [f64; 12] = [1.02, 1.88, 2.34, 2.98, 3.27, 3.26, 3.17, 2.77, 1.92, 1.23, 1.00, 0.61];

fn model(p: &DVector<f64>) -> DVector<f64> {
    let (ka, ke, v) = (p[0], p[1], p[2]);
    DVector::from_iterator(12, TT.iter().map(|&t|
        D * ka / (v * (ka - ke)) * ((-ke * t).exp() - (-ka * t).exp())))
}

fn jac(p: &DVector<f64>, h: f64) -> DMatrix<f64> {   // forward differences
    let m0 = model(p);
    let mut j = DMatrix::zeros(12, 3);
    for k in 0..3 {
        let mut pe = p.clone();
        let e = h * p[k];
        pe[k] += e;
        j.set_column(k, &((model(&pe) - &m0) / e));
    }
    j
}

// least squares min ||A x - b|| by Householder QR: R x = Q^T b
fn lstsq_qr(a: &DMatrix<f64>, b: &DVector<f64>) -> DVector<f64> {
    let n = a.ncols();
    let qr = a.clone().qr();
    let mut qtb = b.clone();
    qr.q_tr_mul(&mut qtb);
    qr.r().solve_upper_triangular(&qtb.rows(0, n).into_owned()).unwrap()
}

fn rss(p: &DVector<f64>) -> f64 {
    (DVector::from_row_slice(&Y) - model(p)).norm_squared()
}

fn main() {
    let y = DVector::from_row_slice(&Y);
    let mut p = DVector::from_vec(vec![0.5, 0.1, 10.0]);
    for k in 1..=30 {
        let r = &y - model(&p);
        let delta = lstsq_qr(&jac(&p, 1e-7), &r);   // linear least squares by QR
        p += &delta;
        println!("{:2} {:10.6} {:10.6} {:10.4} {:12.8}", k, p[0], p[1], p[2], rss(&p));
        if delta.norm() < 1e-8 * p.norm() { break; }
    }
    let j = jac(&p, 1e-7);
    let cond = |m: &DMatrix<f64>| { let s = m.singular_values(); s.max() / s.min() };
    println!("{} {}", cond(&j), cond(&(j.transpose() * &j)));
}
