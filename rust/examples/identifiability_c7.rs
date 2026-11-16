// C7: fixes: compare sampling designs; a prior on log ka as an extra residual.
#![allow(dead_code)]
use nalgebra::{DMatrix, DVector};

// One-compartment model with oral absorption, dose D = 100 mg (mg/L)
const D: f64 = 100.0; // oral dose, mg
fn bateman(t: f64, ka: f64, ke: f64, v: f64) -> f64 {
    D * ka / (v * (ka - ke)) * ((-ke * t).exp() - (-ka * t).exp())
}
const TT: [f64; 12] = [0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 12.0];
const Y: [f64; 12] = [1.02, 1.88, 2.34, 2.98, 3.27, 3.26, 3.17, 2.77, 1.92, 1.23, 1.00, 0.61];

// Central-difference Jacobian of f at p (like numDeriv::jacobian, less refined)
fn jacobian(f: &dyn Fn(&[f64]) -> Vec<f64>, p: &[f64]) -> DMatrix<f64> {
    let n = f(p).len();
    let mut j = DMatrix::zeros(n, p.len());
    for k in 0..p.len() {
        let h = 1e-5 * p[k].abs().max(1.0);
        let (mut pp, mut pm) = (p.to_vec(), p.to_vec());
        pp[k] += h;
        pm[k] -= h;
        let (fp, fm) = (f(&pp), f(&pm));
        for i in 0..n { j[(i, k)] = (fp[i] - fm[i]) / (2.0 * h); }
    }
    j
}

// Levenberg-Marquardt for min sum r(p)^2, hand-written (R: nls / nls.lm, Python: least_squares)
fn levmar(r: &dyn Fn(&[f64]) -> Vec<f64>, p0: &[f64]) -> Vec<f64> {
    let ss = |p: &[f64]| -> f64 { let v = r(p); if v.iter().all(|x| x.is_finite()) { v.iter().map(|x| x * x).sum() } else { f64::INFINITY } };
    let (mut p, mut mu) = (p0.to_vec(), 1e-3);
    for _ in 0..500 {
        let j = jacobian(r, &p);                     // Jacobian of the residuals
        let rv = DVector::from_vec(r(&p));
        let jtj = j.transpose() * &j;
        let g = j.transpose() * &rv;
        loop {
            let a = &jtj + DMatrix::from_diagonal(&jtj.diagonal()) * mu;
            let step = a.lu().solve(&(-&g)).unwrap();
            let pn: Vec<f64> = p.iter().zip(step.iter()).map(|(a, b)| a + b).collect();
            if ss(&pn) <= ss(&p) {
                let small = step.norm() < 1e-12 * (1.0 + DVector::from_vec(p.clone()).norm());
                p = pn;
                mu = (mu / 10.0).max(1e-12);
                if small { return p; }
                break;
            }
            mu *= 10.0;
            if mu > 1e12 { return p; }
        }
    }
    p
}

// Residuals y - C(t) for the one-compartment model at the sample times ts
fn resid(ts: &[f64], ys: &[f64]) -> impl Fn(&[f64]) -> Vec<f64> {
    let (ts, ys) = (ts.to_vec(), ys.to_vec());
    move |p: &[f64]| ts.iter().zip(&ys).map(|(ti, yi)| yi - bateman(*ti, p[0], p[1], p[2])).collect()
}
fn sse(r: &dyn Fn(&[f64]) -> Vec<f64>, p: &[f64]) -> f64 { r(p).iter().map(|x| x * x).sum() }

// predicted SE of log ka, log ke, log V for sample times ts, at a guessed theta
fn se_design(ts: &[f64], p: [f64; 3], sigma: f64) -> Vec<f64> {
    let f = |lp: &[f64]| -> Vec<f64> { ts.iter().map(|t| bateman(*t, lp[0].exp(), lp[1].exp(), lp[2].exp())).collect() };
    let lp: Vec<f64> = p.iter().map(|x| x.ln()).collect();
    let j = jacobian(&f, &lp);
    let c = (j.transpose() * &j).try_inverse().unwrap();
    (0..3).map(|i| sigma * c[(i, i)].sqrt()).collect()
}

fn main() {
    let designs: [(&str, [f64; 6]); 3] = [("late", [3.0, 4.0, 6.0, 8.0, 10.0, 12.0]),
        ("early", [0.5, 1.0, 2.0, 4.0, 8.0, 12.0]), ("peak", [1.0, 2.0, 3.0, 4.0, 6.0, 12.0])];
    for (name, ts) in designs.iter() {
        let se = se_design(ts, [1.0, 0.2, 20.0], 0.1);
        println!("{:6} {:.3} {:.3} {:.3}", name, se[0], se[1], se[2]);
    }
    // (b) prior on log ka (centre log 1.2, SD tau), sparse data, noise SD 0.1
    for tau in [0.1, 0.25, 0.5, 1.0, 2.0] {
        let res = |lp: &[f64]| -> Vec<f64> {
            let mut r: Vec<f64> = TT[6..].iter().zip(&Y[6..])
                .map(|(t, y)| (y - bateman(*t, lp[0].exp(), lp[1].exp(), lp[2].exp())) / 0.1).collect();
            r.push((lp[0] - 1.2_f64.ln()) / tau);
            r
        };
        let lp = levmar(&res, &[0.0, 0.2_f64.ln(), 20.0_f64.ln()]);
        let ja = jacobian(&res, &lp);
        let c = (ja.transpose() * &ja).try_inverse().unwrap();
        println!("tau {:.2}  ka {:.3}  se_logka {:.3}", tau, lp[0].exp(), c[(0, 0)].sqrt());
    }
}
