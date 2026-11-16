// C3: structural non-identifiability: F and V together; the flip-flop in a multistart.
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
    // (a) F and V both free: C depends on V / F only
    let r4 = |p: &[f64]| -> Vec<f64> { T.iter().zip(&Y).map(|(t, y)| y - pk(*t, p[0], p[1], p[2] / p[3])).collect() };
    let p4 = [1.0151_f64, 0.1906, 20.22, 1.0];
    let lp4: Vec<f64> = p4.iter().map(|x| x.ln()).collect();
    let flog = |lp: &[f64]| -> Vec<f64> { T.iter().map(|t| pk(*t, lp[0].exp(), lp[1].exp(), (lp[2] - lp[3]).exp())).collect() };
    let svd = jacobian(&flog, &lp4).svd(false, true);
    let mut idx: Vec<usize> = (0..4).collect();          // nalgebra does not promise an order
    idx.sort_by(|a, b| svd.singular_values[*b].partial_cmp(&svd.singular_values[*a]).unwrap());
    let d: Vec<String> = idx.iter().map(|&i| format!("{:.3e}", svd.singular_values[i])).collect();
    println!("singular values: {}", d.join(" "));
    let vt = svd.v_t.unwrap();
    let v4 = vt.row(idx[3]) * vt[(idx[3], 2)].signum();
    println!("v4: {:.3} {:.3} {:.3} {:.3}", v4[0], v4[1], v4[2], v4[3]);
    let p = levmar(&r4, &[0.5, 0.1, 10.0, 0.8]);
    println!("LM with F free: ka {:.4} ke {:.4} V {:.4} F {:.4}  V/F {:.4}", p[0], p[1], p[2], p[3], p[2] / p[3]);
    // (b) multistart
    let r = resid(&T, &Y);
    for (ka0, ke0) in [(0.5, 0.1), (2.0, 0.1), (0.1, 0.5), (2.0, 0.5)] {
        let th = levmar(&r, &[ka0, ke0, 10.0]);
        println!("start ka={:.1} ke={:.1} -> ka={:.4} ke={:.4} V={:.3} SSE={:.6}",
                 ka0, ke0, th[0], th[1], th[2], sse(&r, &th));
    }
}
