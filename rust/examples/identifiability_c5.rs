// C5: profile likelihood of ka (others re-fitted), and profile vs Wald intervals.
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

// fix ka, re-fit ke and V; return the minimal SSE
fn profile_ka(ka: f64, ts: &[f64], ys: &[f64]) -> f64 {
    let r = |q: &[f64]| -> Vec<f64> { ts.iter().zip(ys).map(|(t, y)| y - bateman(*t, ka, q[0], q[1])).collect() };
    let q = levmar(&r, &[0.2, 20.0]);
    sse(&r, &q)
}

// bisection for a sign change of f on [a, b] (R: uniroot, Python: brentq)
fn bisect(f: &dyn Fn(f64) -> f64, mut a: f64, mut b: f64) -> f64 {
    let fa = f(a);
    for _ in 0..60 {
        let m = 0.5 * (a + b);
        if (f(m) > 0.0) == (fa > 0.0) { a = m } else { b = m }
    }
    0.5 * (a + b)
}

fn main() {
    let grid = [0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 5.0, 10.0, 30.0];
    for start in [0, 6] {                      // all samples; samples from 3 h
        let (ts, ys) = (&TT[start..], &Y[start..]);
        let r = resid(ts, ys);
        let th = levmar(&r, &[1.0, 0.2, 20.0]);
        let (smin, s2) = (sse(&r, &th), sse(&r, &th) / (ts.len() - 3) as f64);
        let row: Vec<String> = grid.iter().map(|&k| format!("{:.2}", (profile_ka(k, ts, ys) - smin) / s2)).collect();
        println!("delta: {}", row.join(" "));
    }
    // (b) profile interval vs Wald intervals, all samples
    let r = resid(&TT, &Y);
    let th = levmar(&r, &[1.0, 0.2, 20.0]);
    let (smin, s2) = (sse(&r, &th), sse(&r, &th) / 9.0);
    let crit = |k: f64| (profile_ka(k, &TT, &Y) - smin) / s2 - 3.84;
    let (lo, hi) = (bisect(&crit, 0.3, th[0]), bisect(&crit, th[0], 30.0));
    let j = jacobian(&r, &th);
    let se = ((j.transpose() * &j).try_inverse().unwrap()[(0, 0)] * s2).sqrt();
    println!("profile  {:.4} {:.4}", lo, hi);
    println!("wald     {:.4} {:.4}", th[0] - 1.96 * se, th[0] + 1.96 * se);
    println!("wald_log {:.4} {:.4}", th[0] * (-1.96 * se / th[0]).exp(), th[0] * (1.96 * se / th[0]).exp());
}
