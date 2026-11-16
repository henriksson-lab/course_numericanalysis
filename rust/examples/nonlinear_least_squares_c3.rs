// C3: Gauss-Newton overshoots from a bad start; Levenberg-Marquardt gets there
use nalgebra::{DMatrix, DVector};

const D: f64 = 100.0;
const TT: [f64; 12] = [0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 12.0];
const Y: [f64; 12] = [1.02, 1.88, 2.34, 2.98, 3.27, 3.26, 3.17, 2.77, 1.92, 1.23, 1.00, 0.61];

fn model(p: &DVector<f64>) -> DVector<f64> {
    let (ka, ke, v) = (p[0], p[1], p[2]);
    DVector::from_iterator(12, TT.iter().map(|&t|
        D * ka / (v * (ka - ke)) * ((-ke * t).exp() - (-ka * t).exp())))
}

fn jac(p: &DVector<f64>, h: f64) -> DMatrix<f64> {
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

// scientific notation as R prints it: 1e-01 rather than Rust's 1e-1
fn sci(x: f64) -> String {
    let s = format!("{:.0e}", x);
    let (m, e) = s.split_once('e').unwrap();
    let e: i32 = e.parse().unwrap();
    format!("{}e{}{:02}", m, if e < 0 { '-' } else { '+' }, e.abs())
}

fn lev_mar(mut p: DVector<f64>, mut lambda: f64, maxit: usize, tol: f64) -> DVector<f64> {
    let y = DVector::from_row_slice(&Y);
    let mut rss_old = rss(&p);
    for k in 1..=maxit {
        let j = jac(&p, 1e-7);
        let r = &y - model(&p);
        // augmented least-squares problem [J; sqrt(lambda) Dg] delta = [r; 0]
        let mut a = DMatrix::zeros(15, 3);
        a.rows_mut(0, 12).copy_from(&j);
        for c in 0..3 {
            a[(12 + c, c)] = lambda.sqrt() * j.column(c).norm();   // Marquardt's scaling
        }
        let mut b = DVector::zeros(15);
        b.rows_mut(0, 12).copy_from(&r);
        let delta = lstsq_qr(&a, &b);
        let p_new = &p + &delta;
        let rss_new = rss(&p_new);
        if rss_new.is_finite() && rss_new < rss_old {       // accept
            p = p_new; rss_old = rss_new; lambda /= 10.0;
            println!("{:2} {:10.6} {:10.6} {:10.4} {:12.8} {:>8}", k, p[0], p[1], p[2], rss_old, sci(lambda));
            if delta.norm() < tol * p.norm() { break; }
        } else {                                            // reject
            lambda *= 10.0;
        }
    }
    p
}

fn main() {
    let y = DVector::from_row_slice(&Y);
    let p0 = DVector::from_vec(vec![3.0, 0.05, 50.0]);
    let p1 = &p0 + lstsq_qr(&jac(&p0, 1e-7), &(&y - model(&p0)));   // one Gauss-Newton step
    println!("{:.8} {:.8} {:.8}  RSS {:e}", p1[0], p1[1], p1[2], rss(&p1));
    let p = lev_mar(p0, 0.01, 200, 1e-8);
    println!("{:.7} {:.7} {:.7}", p[0], p[1], p[2]);
}
