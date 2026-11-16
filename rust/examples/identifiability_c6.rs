// C6: a sloppy model: eigenvalues of J^T J for the two-compartment PK model.
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

// Two-compartment model with oral absorption: closed-form tri-exponential solution
fn pk2(t: f64, p: &[f64]) -> f64 {
    let (ka, k10, k12, k21, v) = (p[0], p[1], p[2], p[3], p[4]);
    let s = k10 + k12 + k21;
    let q = (s * s - 4.0 * k21 * k10).sqrt();
    let (al, be) = ((s + q) / 2.0, (s - q) / 2.0);       // the two disposition rates
    ka * 100.0 / v * ((k21 - al) / ((ka - al) * (be - al)) * (-al * t).exp()
        + (k21 - be) / ((ka - be) * (al - be)) * (-be * t).exp()
        + (k21 - ka) / ((al - ka) * (be - ka)) * (-ka * t).exp())
}

fn sloppy(ts: &[f64]) {
    let lp2: Vec<f64> = [1.0_f64, 0.2, 0.3, 0.15, 20.0].iter().map(|x| x.ln()).collect();
    let f = |lp: &[f64]| -> Vec<f64> {
        let p: Vec<f64> = lp.iter().map(|x| x.exp()).collect();
        ts.iter().map(|t| pk2(*t, &p)).collect()
    };
    let j = jacobian(&f, &lp2);
    let eig = (j.transpose() * &j).symmetric_eigen();
    let mut idx: Vec<usize> = (0..5).collect();          // sort largest first
    idx.sort_by(|a, b| eig.eigenvalues[*b].partial_cmp(&eig.eigenvalues[*a]).unwrap());
    let lam: Vec<f64> = idx.iter().map(|&i| eig.eigenvalues[i]).collect();
    let fmt = |v: Vec<String>| v.join(" ");
    println!("eigenvalues: {}", fmt(lam.iter().map(|x| format!("{:.3e}", x)).collect()));
    println!("decades: {}", fmt(lam.iter().map(|x| format!("{:.2}", (lam[0] / x).log10())).collect()));
    println!("SE along eigenvectors (sigma = 0.1): {}", fmt(lam.iter().map(|x| format!("{:.3}", 0.1 / x.sqrt())).collect()));
    let v5 = eig.eigenvectors.column(idx[4]) * eig.eigenvectors[(0, idx[4])].signum();
    println!("sloppiest (ka, k10, k12, k21, V): {}", fmt(v5.iter().map(|x| format!("{:.2}", x)).collect()));
}

fn main() {
    let p = [1.0, 0.2, 0.3, 0.15, 20.0];
    println!("{:.4} {:.4} {:.4}", pk2(1.0, &p), pk2(6.0, &p), pk2(24.0, &p));
    sloppy(&TT);
    let mut longer = TT.to_vec();
    longer.extend([24.0, 36.0, 48.0]);
    sloppy(&longer);
}
