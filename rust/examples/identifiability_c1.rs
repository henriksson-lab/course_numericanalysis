// C1: standard errors three ways, and a check of B3.
#![allow(dead_code)]
use nalgebra::{DMatrix, DVector};

// One-compartment model with oral absorption, dose D = 100 mg (mg/L)
fn pk(t: f64, ka: f64, ke: f64, v: f64) -> f64 {
    100.0 * ka / (v * (ka - ke)) * ((-ke * t).exp() - (-ka * t).exp())
}
const T: [f64; 12] = [0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 12.0];
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

// Residuals y - C(t) for the one-compartment model at the sample times tt
fn resid(tt: &[f64], yy: &[f64]) -> impl Fn(&[f64]) -> Vec<f64> {
    let (tt, yy) = (tt.to_vec(), yy.to_vec());
    move |p: &[f64]| tt.iter().zip(&yy).map(|(ti, yi)| yi - pk(*ti, p[0], p[1], p[2])).collect()
}
fn sse(r: &dyn Fn(&[f64]) -> Vec<f64>, p: &[f64]) -> f64 { r(p).iter().map(|x| x * x).sum() }

fn main() {
    let r = resid(&T, &Y);
    let th = levmar(&r, &[0.5, 0.1, 10.0]);
    println!("estimate: {:.5} {:.5} {:.5}", th[0], th[1], th[2]);
    let s2 = sse(&r, &th) / (T.len() - 3) as f64;
    let j = jacobian(&r, &th);                        // = -J; J^T J is the same
    let c_solve = (j.transpose() * &j).try_inverse().unwrap() * s2;
    let rr = j.clone().qr().r();                      // (J^T J)^-1 = R^-1 R^-T
    let rinv = rr.try_inverse().unwrap();
    let c_qr = &rinv * rinv.transpose() * s2;
    for (name, c) in [("se_solve", &c_solve), ("se_qr   ", &c_qr)] {
        println!("{} {:.5} {:.5} {:.5}", name, c[(0, 0)].sqrt(), c[(1, 1)].sqrt(), c[(2, 2)].sqrt());
    }
    println!("correlation:");
    for i in 0..3 {
        let row: Vec<String> = (0..3)
            .map(|k| format!("{:7.3}", c_qr[(i, k)] / (c_qr[(i, i)] * c_qr[(k, k)]).sqrt())).collect();
        println!("{}", row.join(" "));
    }
    // (d) check B3: C = A exp(-k t), A = 10, k = 0.5, t = 0, 1, 2, sigma = 0.2
    let (a, k) = (10.0_f64, 0.5_f64);
    let jb = DMatrix::from_fn(3, 2, |i, c| {
        let t = i as f64;
        if c == 0 { (-k * t).exp() } else { -a * t * (-k * t).exp() }
    });
    let cb = (jb.transpose() * &jb).try_inverse().unwrap() * 0.04;
    println!("B3: SE {:.4} {:.4}, correlation {:.4}", cb[(0, 0)].sqrt(), cb[(1, 1)].sqrt(),
             cb[(0, 1)] / (cb[(0, 0)] * cb[(1, 1)]).sqrt());
}
