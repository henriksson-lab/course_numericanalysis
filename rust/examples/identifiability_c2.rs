// C2: the Hessian of the SSE by central differences; SEs from it.
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

// Central-difference Hessian of a scalar function (R: numDeriv::hessian)
fn hessian(f: &dyn Fn(&[f64]) -> f64, p: &[f64]) -> DMatrix<f64> {
    let n = p.len();
    let h: Vec<f64> = p.iter().map(|x| 1e-4 * x.abs().max(1.0)).collect();
    let fe = |di: usize, si: f64, dk: usize, sk: f64| {
        let mut q = p.to_vec();
        q[di] += si * h[di];
        q[dk] += sk * h[dk];
        f(&q)
    };
    DMatrix::from_fn(n, n, |i, k| {
        (fe(i, 1.0, k, 1.0) - fe(i, 1.0, k, -1.0) - fe(i, -1.0, k, 1.0) + fe(i, -1.0, k, -1.0))
            / (4.0 * h[i] * h[k])
    })
}

fn main() {
    let r = resid(&T, &Y);
    let th = levmar(&r, &[0.5, 0.1, 10.0]);
    let s2 = sse(&r, &th) / (T.len() - 3) as f64;
    let j = jacobian(&r, &th);
    let f = |p: &[f64]| sse(&r, p);
    let h = hessian(&f, &th);
    let gn = (j.transpose() * &j) * 2.0;
    let rel = (0..9).map(|i| ((h[i] - gn[i]) / h[i]).abs()).fold(0.0, f64::max);
    println!("largest relative difference H vs 2 J^T J: {:.4}", rel);
    let c = h.try_inverse().unwrap() * (2.0 * s2);
    println!("SE from H: {:.6} {:.6} {:.6}", c[(0, 0)].sqrt(), c[(1, 1)].sqrt(), c[(2, 2)].sqrt());
    // (c) known sigma = 0.1: minimise nll = SSE / (2 sigma^2) (LM on r / sigma), invert its Hessian
    let nll = |p: &[f64]| sse(&r, p) / (2.0 * 0.01);
    for start in [[1.0, 0.2, 20.0], [0.5, 0.1, 10.0]] {
        let p = levmar(&r, &start);
        let ci = hessian(&nll, &p).try_inverse().unwrap();
        println!("{:.4} {:.4} {:.4}  SE {:.4} {:.4} {:.4}", p[0], p[1], p[2],
                 ci[(0, 0)].sqrt(), ci[(1, 1)].sqrt(), ci[(2, 2)].sqrt());
    }
}
