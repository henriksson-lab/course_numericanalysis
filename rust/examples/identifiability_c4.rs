// C4: the SVD of the Jacobian in log-parameters, both designs.
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

fn svd_report(ts: &[f64], ys: &[f64]) {
    let r = resid(ts, ys);
    let th = levmar(&r, &[1.0, 0.2, 20.0]);
    let s = (sse(&r, &th) / (ts.len() - 3) as f64).sqrt();
    let model = |p: &[f64]| -> Vec<f64> { ts.iter().map(|t| bateman(*t, p[0], p[1], p[2])).collect() };
    let j = jacobian(&model, &th);
    let jl = DMatrix::from_fn(ts.len(), 3, |i, k| j[(i, k)] * th[k]);  // d C / d log(theta)
    let kappa = |m: &DMatrix<f64>| { let d = m.singular_values(); d.max() / d.min() };
    let svd = jl.clone().svd(false, true);
    let mut idx: Vec<usize> = (0..3).collect();
    idx.sort_by(|a, b| svd.singular_values[*b].partial_cmp(&svd.singular_values[*a]).unwrap());
    let d: Vec<f64> = idx.iter().map(|&i| svd.singular_values[i]).collect();
    println!("singular values: {:.4} {:.4} {:.4}", d[0], d[1], d[2]);
    println!("kappa(J) = {:.1}  kappa(J_log) = {:.3}", kappa(&j), kappa(&jl));
    println!("SE along v_i (log units): {:.3} {:.3} {:.3}", s / d[0], s / d[1], s / d[2]);
    let vt = svd.v_t.unwrap();
    let v3 = vt.row(idx[2]) * vt[(idx[2], 0)].signum();
    println!("v3 (log ka, log ke, log V): {:.3} {:.3} {:.3}", v3[0], v3[1], v3[2]);
}

fn main() {
    svd_report(&TT, &Y);
    svd_report(&TT[6..], &Y[6..]);          // samples from 3 h
}
