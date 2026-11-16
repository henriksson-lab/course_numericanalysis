// C7: multistart: Levenberg-Marquardt (log parameters) from a 6 x 6 grid of starts
use nalgebra::{DMatrix, DVector};

const D: f64 = 100.0;
const TT: [f64; 12] = [0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 12.0];
const Y: [f64; 12] = [1.02, 1.88, 2.34, 2.98, 3.27, 3.26, 3.17, 2.77, 1.92, 1.23, 1.00, 0.61];

fn f_nat(p: &DVector<f64>) -> DVector<f64> {               // p = (ka, ke, V)
    let (ka, ke, v) = (p[0], p[1], p[2]);
    DVector::from_iterator(12, TT.iter().map(|&t|
        D * ka / (v * (ka - ke)) * ((-ke * t).exp() - (-ka * t).exp())))
}
fn f_log(q: &DVector<f64>) -> DVector<f64> { f_nat(&q.map(f64::exp)) }   // q = log(p)

fn jac(f: &dyn Fn(&DVector<f64>) -> DVector<f64>, p: &DVector<f64>) -> DMatrix<f64> {
    let m0 = f(p);
    let mut j = DMatrix::zeros(m0.len(), p.len());
    for k in 0..p.len() {
        let mut pe = p.clone();
        let e = 1e-7 * p[k].abs().max(1.0);
        pe[k] += e;
        j.set_column(k, &((f(&pe) - &m0) / e));
    }
    j
}

fn lstsq_qr(a: &DMatrix<f64>, b: &DVector<f64>) -> Option<DVector<f64>> {
    let n = a.ncols();
    let qr = a.clone().qr();
    let mut qtb = b.clone();
    qr.q_tr_mul(&mut qtb);
    qr.r().solve_upper_triangular(&qtb.rows(0, n).into_owned())   // None if R is singular
}

fn rss(f: &dyn Fn(&DVector<f64>) -> DVector<f64>, p: &DVector<f64>) -> f64 {
    (DVector::from_row_slice(&Y) - f(p)).norm_squared()
}

// returns (parameters, RSS, iterations)
fn lev_mar(f: &dyn Fn(&DVector<f64>) -> DVector<f64>, mut p: DVector<f64>)
    -> (DVector<f64>, f64, usize) {
    let (n, m) = (p.len(), Y.len());
    let y = DVector::from_row_slice(&Y);
    let (mut lambda, mut rss_old, mut k) = (0.01f64, rss(f, &p), 0);
    while k < 200 {
        k += 1;
        let j = jac(f, &p);
        let mut a = DMatrix::zeros(m + n, n);
        a.rows_mut(0, m).copy_from(&j);
        for c in 0..n { a[(m + c, c)] = lambda.sqrt() * j.column(c).norm(); }
        let mut b = DVector::zeros(m + n);
        b.rows_mut(0, m).copy_from(&(&y - f(&p)));
        let delta = match lstsq_qr(&a, &b) { Some(d) => d, None => break };
        if !delta.iter().all(|v| v.is_finite()) { break; }
        if delta.norm() < 1e-8 * (1.0 + p.norm()) { break; }     // converged
        let p_new = &p + &delta;
        let rss_new = rss(f, &p_new);
        if rss_new.is_finite() && rss_new < rss_old {
            p = p_new; rss_old = rss_new; lambda /= 10.0;
        } else {
            lambda *= 10.0;
        }
    }
    (p, rss_old, k)
}

// like C's %.4g
fn g4(x: f64) -> String {
    if !x.is_finite() || x == 0.0 { return format!("{}", x); }
    let e = x.abs().log10().floor() as i32;
    if e < -4 || e >= 4 { format!("{:.3e}", x) } else {
        let s = format!("{:.*}", (3 - e).max(0) as usize, x);
        if s.contains('.') { s.trim_end_matches('0').trim_end_matches('.').to_string() } else { s }
    }
}

fn main() {
    let g: [f64; 6] = [0.03, 0.1, 0.3, 1.0, 3.0, 10.0];
    let mut ends: Vec<(String, usize)> = Vec::new();
    for &ke0 in g.iter() {
        for &ka0 in g.iter() {                    // 36 starts, ke0 = 1.5 g: never ka = ke
            let q0 = DVector::from_vec(vec![ka0.ln(), (1.5 * ke0).ln(), 10.0f64.ln()]);
            let (q, rss, _) = lev_mar(&f_log, q0);
            let p = q.map(f64::exp);
            let key = format!("{}  {}  {}  {}", g4(p[0]), g4(p[1]), g4(p[2]), g4(rss));
            match ends.iter_mut().find(|(e, _)| *e == key) {
                Some((_, n)) => *n += 1,
                None => ends.push((key, 1)),
            }
        }
    }
    ends.sort_by(|a, b| b.1.cmp(&a.1));
    for (end, n) in ends { println!("{:2}   {}", n, end); }
}
