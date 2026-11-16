// C6: interpolation vs least-squares fitting of noisy readings.
use nalgebra::{DMatrix, DVector};

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
fn max_abs_diff(p: &[f64], q: &[f64]) -> f64 {
    p.iter().zip(q).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max)
}


// least-squares polynomial of degree d in s = (t - 6) / 6 (rescaled, like poly())
fn ls_fit(t: &[f64], y: &[f64], d: usize) -> DVector<f64> {
    let a = DMatrix::from_fn(t.len(), d + 1, |i, j| ((t[i] - 6.0) / 6.0).powi(j as i32));
    a.svd(true, true).solve(&DVector::from_column_slice(y), 1e-12).unwrap()
}
fn ls_eval(cf: &DVector<f64>, t: f64) -> f64 {
    let s = (t - 6.0) / 6.0;
    cf.iter().rev().fold(0.0, |acc, &c| acc * s + c)
}

fn main() {
    let t12: Vec<f64> = (0..=12).map(|i| i as f64).collect();
    let y12 = [0.0375, 0.1095, 0.1959, 0.4207, 0.6259, 0.8279, 1.0188, 1.1207,
               1.1673, 1.1736, 1.2210, 1.2037, 1.1857]; // hourly, noise sd 0.02
    let tt: Vec<f64> = (0..=1200).map(|i| i as f64 * 0.01).collect();
    let truth: Vec<f64> = tt.iter().map(|&t| od(t)).collect();
    let p_int = interp(&t12, &y12, &tt); // (a) degree 12
    let (lo, hi) = p_int.iter().fold((f64::INFINITY, f64::NEG_INFINITY),
                                     |(lo, hi), &v| (lo.min(v), hi.max(v)));
    println!("interpolant: range {lo:.3} to {hi:.3}, max error {:.4}",
             max_abs_diff(&p_int, &truth));
    for d in 2..=8 { // (b) least squares
        let cf = ls_fit(&t12, &y12, d);
        let rss: f64 = t12.iter().zip(&y12).map(|(&t, &y)| (y - ls_eval(&cf, t)).powi(2)).sum();
        let p_fit: Vec<f64> = tt.iter().map(|&t| ls_eval(&cf, t)).collect();
        println!("degree {d}: residual sd {:.4}, max error {:.4}",
                 (rss / (13 - d - 1) as f64).sqrt(), max_abs_diff(&p_fit, &truth));
    }
    // (c) piecewise linear
    let p_lin: Vec<f64> = tt.iter().map(|&z| {
        let i = (1..13).find(|&i| z <= t12[i]).unwrap();
        y12[i - 1] + (y12[i] - y12[i - 1]) * (z - t12[i - 1])
    }).collect();
    println!("{}", max_abs_diff(&p_lin, &truth));
    // (d) check without truth
    let cf5 = ls_fit(&t12, &y12, 5);
    let p5: Vec<f64> = tt.iter().map(|&t| ls_eval(&cf5, t)).collect();
    println!("{}", max_abs_diff(&p_int, &p5));
}
