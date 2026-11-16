use nalgebra::{DMatrix, DVector};
fn pk(t: f64, p: &[f64]) -> f64 { let (ka, ke, v) = (p[0], p[1], p[2]);
    100.0 * ka / (v * (ka - ke)) * ((-ke * t).exp() - (-ka * t).exp()) }
fn jac(t: &[f64], p: &[f64]) -> DMatrix<f64> {
    let (ka, ke, v) = (p[0], p[1], p[2]);
    DMatrix::from_fn(t.len(), 3, |i, j| {
        let a = 100.0 * ka / (v * (ka - ke));
        let e = (-ke * t[i]).exp() - (-ka * t[i]).exp();
        match j { 0 => -100.0 * ke / (v * (ka - ke).powi(2)) * e + a * t[i] * (-ka * t[i]).exp(),
                  1 => 100.0 * ka / (v * (ka - ke).powi(2)) * e - a * t[i] * (-ke * t[i]).exp(),
                  _ => -a * e / v } })
}
fn lm(t: &[f64], y: &[f64], p0: &[f64]) -> Vec<f64> {
    let mut p = p0.to_vec(); let mut mu = 1e-3;
    let sse = |p: &[f64]| t.iter().zip(y).map(|(ti, yi)| (yi - pk(*ti, p)).powi(2)).sum::<f64>();
    for _ in 0..200 {
        let j = jac(t, &p);
        let r = DVector::from_fn(t.len(), |i, _| y[i] - pk(t[i], &p));
        let jtj = j.transpose() * &j; let g = j.transpose() * &r;
        loop {
            let a = &jtj + DMatrix::from_diagonal(&jtj.diagonal()) * mu;
            let dp = a.cholesky().unwrap().solve(&g);
            let pn: Vec<f64> = p.iter().zip(dp.iter()).map(|(a, b)| a + b).collect();
            if sse(&pn) < sse(&p) { mu /= 10.0; let done = dp.norm() < 1e-10 * (1.0 + DVector::from_vec(p.clone()).norm()); p = pn; if done {return p;} break; }
            mu *= 10.0; if mu > 1e10 { return p; }
        }
    }
    p
}
fn main() {
    let t = [0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 24.0];
    let y = [0.98, 1.89, 2.69, 3.47, 3.39, 3.00, 2.77, 1.98, 1.35, 0.52, 0.28];
    let p = lm(&t, &y, &[0.5, 0.1, 10.0]);
    println!("{:?}", p);
    let j = jac(&t, &p);
    let sse: f64 = t.iter().zip(&y).map(|(ti, yi)| (yi - pk(*ti, &p)).powi(2)).sum();
    let s2 = sse / 8.0;
    let cov = (j.transpose() * &j).try_inverse().unwrap() * s2;
    println!("{:.6}", cov);
}
