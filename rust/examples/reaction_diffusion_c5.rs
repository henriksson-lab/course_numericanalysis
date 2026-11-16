// C5: decay lengths from 20 synthetic Fgf8-like profiles (one per embryo),
// fitted by a straight line through log c; then k = D / lambda^2.

/// Park-Miller generator, the same numbers as in R and Python.
fn lcg(n: usize, seed: u64) -> Vec<f64> {
    let mut s = seed;
    (0..n).map(|_| { s = (16807 * s) % 2147483647; s as f64 / 2147483647.0 }).collect()
}

/// Least-squares slope of y against x.
fn slope(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len() as f64;
    let (mx, my) = (x.iter().sum::<f64>() / n, y.iter().sum::<f64>() / n);
    let sxy: f64 = x.iter().zip(y).map(|(a, b)| (a - mx) * (b - my)).sum();
    let sxx: f64 = x.iter().map(|a| (a - mx) * (a - mx)).sum();
    sxy / sxx
}

// scientific notation as R and Python print it: 1.26e-03 rather than Rust's 1.26e-3
fn sci(x: f64, digits: usize) -> String {
    let s = format!("{:.*e}", digits, x);
    let (m, e) = s.split_once('e').unwrap();
    let e: i32 = e.parse().unwrap();
    format!("{}e{}{:02}", m, if e < 0 { '-' } else { '+' }, e.abs())
}

fn main() {
    let dd = 50.0; // µm^2/s (approximate FCS value)
    let x: Vec<f64> = (0..=20).map(|i| 20.0 * i as f64).collect(); // 21 positions
    let nx = x.len();
    let z: Vec<f64> = lcg(20 * (nx + 1), 2).iter().map(|v| v - 0.5).collect();
    for bg in [0.0, 0.05] {
        let lf: Vec<f64> = (0..20).map(|e| {
            let ze = &z[e * (nx + 1)..(e + 1) * (nx + 1)]; // column e
            let lam_true = 200.0 * (1.0 + 0.2 * ze[0]);
            let logc: Vec<f64> = (0..nx)
                .map(|i| ((-x[i] / lam_true).exp() * (1.0 + 0.2 * ze[i + 1]) + bg).ln())
                .collect();
            -1.0 / slope(&x, &logc)
        }).collect();
        let mean = lf.iter().sum::<f64>() / 20.0;
        let sd = (lf.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / 19.0).sqrt();
        let kmean = lf.iter().map(|v| dd / (v * v)).sum::<f64>() / 20.0;
        let life = lf.iter().map(|v| v * v / dd).sum::<f64>() / 20.0 / 60.0;
        println!("background {:.2}: lambda = {:.1} +- {:.1} µm, k = {} /s, lifetime {:.1} min",
                 bg, mean, sd, sci(kmean, 2), life);
    }
}
