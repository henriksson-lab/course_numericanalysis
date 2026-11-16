// Background: steady synthesis-diffusion-degradation gradient on a grid with
// ghost points, solved by Thomas, checked against the exact cosh formula.

fn thomas(a: &[f64], b: &[f64], c: &[f64], d: &[f64]) -> Vec<f64> {
    let n = d.len();
    let (mut b, mut d) = (b.to_vec(), d.to_vec()); // copies: inputs stay unchanged
    for i in 1..n {
        let w = a[i] / b[i - 1];
        b[i] -= w * c[i - 1];
        d[i] -= w * d[i - 1];
    }
    let mut x = vec![0.0; n];
    x[n - 1] = d[n - 1] / b[n - 1];
    for i in (0..n - 1).rev() {
        x[i] = (d[i] - c[i] * x[i + 1]) / b[i];
    }
    x
}

fn main() {
    let (dd, k, l) = (100.0_f64, 0.01_f64, 500.0_f64); // µm^2/min, 1/min, µm
    let lam = (dd / k).sqrt();
    let j = dd / lam; // decay length 100 µm; c(0) close to 1
    let n = 50;
    let dx = l / n as f64;
    let r = dd / (dx * dx);
    let x: Vec<f64> = (0..=n).map(|i| i as f64 * dx).collect();
    // -D c'' on the grid, ghost points at both ends
    let mut lo = vec![-r; n + 1];
    lo[0] = 0.0;
    lo[n] = -2.0 * r;
    let mut up = vec![-r; n + 1];
    up[0] = -2.0 * r;
    up[n] = 0.0;
    let di: Vec<f64> = vec![2.0 * r + k; n + 1]; // diagonal incl. decay
    let mut src = vec![0.0; n + 1];
    src[0] = 2.0 * j / dx; // flux j enters at x = 0
    let cs = thomas(&lo, &di, &up, &src);
    let err = x.iter().zip(&cs)
        .map(|(xi, ci)| (ci - j * lam / dd * ((l - xi) / lam).cosh() / (l / lam).sinh()).abs())
        .fold(0.0, f64::max);
    println!("{}", err);
    let total: f64 = cs.iter().sum::<f64>() - (cs[0] + cs[n]) / 2.0;
    println!("{}", k * dx * total); // total decay (trapezoid) = influx j
}
