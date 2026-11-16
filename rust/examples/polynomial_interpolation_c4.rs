// C4: max error against the degree n, three sampling designs.
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


fn maxerr(x: &[f64], a: f64, b: f64) -> f64 {
    let tt: Vec<f64> = (0..=4000).map(|i| a + (b - a) * i as f64 / 4000.0).collect();
    let y: Vec<f64> = x.iter().map(|&t| od(t)).collect();
    let truth: Vec<f64> = tt.iter().map(|&t| od(t)).collect();
    max_abs_diff(&interp(x, &y, &tt), &truth)
}
fn equi(n: usize, a: f64, b: f64) -> Vec<f64> {
    (0..=n).map(|i| a + (b - a) * i as f64 / n as f64).collect()
}

fn main() {
    for n in [4usize, 8, 12, 16, 20, 30, 40] {
        let e8 = maxerr(&equi(n, 0.0, 8.0), 0.0, 8.0);
        let e48 = maxerr(&equi(n, 0.0, 48.0), 0.0, 48.0);
        let c48 = maxerr(&cheb(n, 0.0, 48.0), 0.0, 48.0);
        println!("{n:2}  {e8:9.2e}  {e48:9.2e}  {c48:9.2e}");
    }
    // No plotting here: write the table to a file and plot it in R or Python.
}
