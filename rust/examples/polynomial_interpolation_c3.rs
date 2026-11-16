// C3: the Runge phenomenon on the 48 h run, and Runge's function.
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

fn main() {
    // (a) 48 h run, a reading every 4 h
    let t48: Vec<f64> = (0..=12).map(|i| 4.0 * i as f64).collect();
    let y48: Vec<f64> = t48.iter().map(|&t| od(t)).collect();
    let tt: Vec<f64> = (0..=4800).map(|i| i as f64 * 0.01).collect();
    let truth: Vec<f64> = tt.iter().map(|&t| od(t)).collect();
    let pe = interp(&t48, &y48, &tt);
    let imin = (0..pe.len()).min_by(|&i, &j| pe[i].total_cmp(&pe[j])).unwrap();
    println!("equispaced: min {:.4} at {:.2} h, max error {:.4}",
             pe[imin], tt[imin], max_abs_diff(&pe, &truth));
    // (b) 13 Chebyshev times
    let tc = cheb(12, 0.0, 48.0);
    let yc: Vec<f64> = tc.iter().map(|&t| od(t)).collect();
    let pc = interp(&tc, &yc, &tt);
    let pcmin = pc.iter().cloned().fold(f64::INFINITY, f64::min);
    println!("Chebyshev:  min {:.4}, max error {:.4}", pcmin, max_abs_diff(&pc, &truth));
    // (c) Runge's function
    let runge = |x: f64| 1.0 / (1.0 + 25.0 * x * x);
    let xx: Vec<f64> = (0..=2000).map(|i| -1.0 + i as f64 * 0.001).collect();
    let rx: Vec<f64> = xx.iter().map(|&x| runge(x)).collect();
    for n in [5usize, 10, 20] {
        let xe: Vec<f64> = (0..=n).map(|i| -1.0 + 2.0 * i as f64 / n as f64).collect();
        let xc = cheb(n, -1.0, 1.0);
        let ye: Vec<f64> = xe.iter().map(|&x| runge(x)).collect();
        let yc: Vec<f64> = xc.iter().map(|&x| runge(x)).collect();
        println!("n = {n:2}: equispaced {:.4}   Chebyshev {:.4}",
                 max_abs_diff(&interp(&xe, &ye, &xx), &rx),
                 max_abs_diff(&interp(&xc, &yc, &xx), &rx));
    }
}
