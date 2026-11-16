// C6: how the tolerance of the inner ODE solver (Dopri5) affects the fit
use nalgebra::{DMatrix, DVector};
use ode_solvers::continuous_output_model::ContinuousOutputModel;
use ode_solvers::{Dopri5, System, Vector2};

const D: f64 = 100.0;
const TT: [f64; 12] = [0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 12.0];
const Y: [f64; 12] = [1.02, 1.88, 2.34, 2.98, 3.27, 3.26, 3.17, 2.77, 1.92, 1.23, 1.00, 0.61];

struct Pk { ka: f64, ke: f64 }
impl System<f64, Vector2<f64>> for Pk {
    fn system(&self, _t: f64, a: &Vector2<f64>, da: &mut Vector2<f64>) {
        da[0] = -self.ka * a[0];                        // gut
        da[1] = self.ka * a[0] - self.ke * a[1];        // central compartment
    }
}

// concentrations at the sampling times: amount in the central compartment / V
fn pk_ode(p: &DVector<f64>, tol: f64) -> DVector<f64> {
    let mut s = Dopri5::new(Pk { ka: p[0], ke: p[1] }, 0.0, 12.0, 12.0,
                            Vector2::new(D, 0.0), tol, tol);
    let mut dense = ContinuousOutputModel::default();
    if s.integrate_with_continuous_output_model(&mut dense).is_err() {
        return DVector::from_element(12, f64::NAN);
    }
    DVector::from_iterator(12, TT.iter().map(|&t| dense.evaluate(t).unwrap()[1] / p[2]))
}

fn bateman(t: f64, ka: f64, ke: f64, v: f64) -> f64 {
    D * ka / (v * (ka - ke)) * ((-ke * t).exp() - (-ka * t).exp())
}

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

fn lstsq_qr(a: &DMatrix<f64>, b: &DVector<f64>) -> DVector<f64> {
    let n = a.ncols();
    let qr = a.clone().qr();
    let mut qtb = b.clone();
    qr.q_tr_mul(&mut qtb);
    qr.r().solve_upper_triangular(&qtb.rows(0, n).into_owned()).unwrap()
}

fn lev_mar(f: &dyn Fn(&DVector<f64>) -> DVector<f64>, mut p: DVector<f64>)
    -> (DVector<f64>, f64, usize) {
    let (n, m) = (p.len(), Y.len());
    let y = DVector::from_row_slice(&Y);
    let rss = |p: &DVector<f64>| (&y - f(p)).norm_squared();
    let (mut lambda, mut rss_old, mut k) = (0.01f64, rss(&p), 0);
    while k < 200 {
        k += 1;
        let j = jac(f, &p);
        let mut a = DMatrix::zeros(m + n, n);
        a.rows_mut(0, m).copy_from(&j);
        for c in 0..n { a[(m + c, c)] = lambda.sqrt() * j.column(c).norm(); }
        let mut b = DVector::zeros(m + n);
        b.rows_mut(0, m).copy_from(&(&y - f(&p)));
        let delta = lstsq_qr(&a, &b);
        if delta.norm() < 1e-8 * (1.0 + p.norm()) { break; }
        let p_new = &p + &delta;
        let rss_new = rss(&p_new);
        if rss_new.is_finite() && rss_new < rss_old {
            p = p_new; rss_old = rss_new; lambda /= 10.0;
        } else {
            lambda *= 10.0;
        }
    }
    (p, rss_old, k)
}

fn main() {
    let y = DVector::from_row_slice(&Y);
    let rss_exact = |p: &DVector<f64>| -> f64 {
        TT.iter().enumerate().map(|(i, &t)| (y[i] - bateman(t, p[0], p[1], p[2])).powi(2)).sum()
    };
    for tol in [1e-2, 1e-3, 1e-4, 1e-6, 1e-8, 1e-10] {
        let f = |q: &DVector<f64>| pk_ode(&q.map(f64::exp), tol);
        let q0 = DVector::from_vec(vec![0.5f64.ln(), 0.1f64.ln(), 10.0f64.ln()]);
        let (q, rss, k) = lev_mar(&f, q0);
        let p = q.map(f64::exp);
        println!("{:.0e}  {:.7} {:.7} {:.5}  RSS {:.8} (true {:.8})  iter {}",
                 tol, p[0], p[1], p[2], rss, rss_exact(&p), k);
    }
    // (b) noise in the objective along a line in ka
    let p_hat = [1.0150721, 0.1906157, 20.2162104];
    for tol in [1e-3, 1e-6, 1e-10] {
        let mut worst = 0.0f64;
        for i in 0..=400 {
            let p = DVector::from_vec(vec![0.9 + 0.2 * i as f64 / 400.0, p_hat[1], p_hat[2]]);
            let rss_ode = (&y - pk_ode(&p, tol)).norm_squared();
            worst = worst.max((rss_ode - rss_exact(&p)).abs());
        }
        println!("{:.0e}  max |RSS_ode - RSS_exact| = {:.1e}", tol, worst);
    }
    let p1 = DVector::from_vec(vec![1.01 * p_hat[0], p_hat[1], p_hat[2]]);
    println!("{}", rss_exact(&p1) - rss_exact(&DVector::from_row_slice(&p_hat)));
}
