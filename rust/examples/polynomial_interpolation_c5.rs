// C5: monomial coefficients (Vandermonde + LU) against the barycentric formula.
use nalgebra::{DMatrix, DVector};
use std::f64::consts::PI;

const K: f64 = 1.2;
const N0: f64 = 0.05;
const R: f64 = 0.8;

fn od(t: f64) -> f64 {
    K / (1.0 + (K / N0 - 1.0) * (-R * t).exp())
}

fn bary_weights(x: &[f64]) -> Vec<f64> {
    (0..x.len())
        .map(|j| 1.0 / (0..x.len()).filter(|&k| k != j).map(|k| x[j] - x[k]).product::<f64>())
        .collect()
}
fn bary_eval(x: &[f64], y: &[f64], w: &[f64], z: f64) -> f64 {
    let (mut num, mut den) = (0.0, 0.0);
    for j in 0..x.len() {
        let d = z - x[j];
        if d == 0.0 {
            return y[j]; // z is a node
        }
        num += w[j] * y[j] / d;
        den += w[j] / d;
    }
    num / den
}
fn interp(x: &[f64], y: &[f64], zs: &[f64]) -> Vec<f64> {
    let w = bary_weights(x);
    zs.iter().map(|&z| bary_eval(x, y, &w, z)).collect()
}
fn cheb(n: usize, a: f64, b: f64) -> Vec<f64> {
    (0..=n)
        .map(|k| (a + b) / 2.0 - (b - a) / 2.0 * ((2 * k + 1) as f64 * PI / (2 * n + 2) as f64).cos())
        .collect()
}
fn max_abs_diff(p: &[f64], q: &[f64]) -> f64 {
    p.iter().zip(q).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max)
}


fn vander(x: &[f64], m: usize) -> DMatrix<f64> {
    DMatrix::from_fn(x.len(), m, |i, j| x[i].powi(j as i32))
}
fn cond(a: &DMatrix<f64>) -> f64 {
    let s = a.clone().svd(false, false).singular_values;
    s.max() / s.min()
}
// solve V c = y by LU (nalgebra never refuses unless exactly singular), evaluate at xe
fn mono_vals(x: &[f64], y: &[f64], xe: &[f64]) -> Vec<f64> {
    let m = x.len();
    let cf = vander(x, m).lu().solve(&DVector::from_column_slice(y)).unwrap();
    (vander(xe, m) * cf).iter().cloned().collect()
}

fn main() {
    // (a) hours, a reading every 4 h
    let t48: Vec<f64> = (0..=12).map(|i| 4.0 * i as f64).collect();
    let y48: Vec<f64> = t48.iter().map(|&t| od(t)).collect();
    println!("{:.6e}", cond(&vander(&t48, 13)));
    let tt: Vec<f64> = (0..=4800).map(|i| i as f64 * 0.01).collect();
    println!("{:e}", max_abs_diff(&mono_vals(&t48, &y48, &tt), &interp(&t48, &y48, &tt)));
    // (b) Chebyshev times: monomials (hours, rescaled) vs barycentric
    let xx: Vec<f64> = (0..=2000).map(|i| 48.0 * i as f64 / 2000.0).collect();
    let truth: Vec<f64> = xx.iter().map(|&t| od(t)).collect();
    let ss: Vec<f64> = xx.iter().map(|&t| (t - 24.0) / 24.0).collect();
    for n in [8usize, 12, 20, 30, 40, 60, 80] {
        let tc = cheb(n, 0.0, 48.0);
        let y: Vec<f64> = tc.iter().map(|&t| od(t)).collect();
        let s: Vec<f64> = tc.iter().map(|&t| (t - 24.0) / 24.0).collect();
        let w: Vec<f64> = (0..=n)   // Chebyshev weights
            .map(|k| (-1f64).powi(k as i32) * ((2 * k + 1) as f64 * PI / (2 * n + 2) as f64).sin())
            .collect();
        let pb: Vec<f64> = xx.iter().map(|&z| bary_eval(&tc, &y, &w, z)).collect();
        println!("n = {n:2}  kappa {:8.2e}  hours {:8.2e}  rescaled {:8.2e}  bary {:8.2e}",
                 cond(&vander(&s, n + 1)),
                 max_abs_diff(&mono_vals(&tc, &y, &xx), &truth),
                 max_abs_diff(&mono_vals(&s, &y, &ss), &truth),
                 max_abs_diff(&pb, &truth));
    }
}
