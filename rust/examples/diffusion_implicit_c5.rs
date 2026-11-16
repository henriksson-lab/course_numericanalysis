// C5: the embryo between a 20 C and a 27 C stream: steady temperature profile
// (one tridiagonal solve), 10-90 % width, observed order, and the effect of the radius.

/// Thomas algorithm (a[0] and c[n-1] unused), as in Structure: symmetric, banded, sparse.
fn thomas(a: &[f64], b: &[f64], c: &[f64], d: &[f64]) -> Vec<f64> {
    let n = d.len();
    let (mut b, mut d) = (b.to_vec(), d.to_vec());
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

const D: f64 = 1.4e5;   // thermal diffusivity of water, um^2/s
const L: f64 = 500.0;   // embryo length, um

fn tau_of(r: f64) -> f64 { r * r / (2.405 * 2.405 * D) }   // exchange time, s

/// Steady temperature on the nodes 0..N; no flux through the tips.
fn embryo_steady(n: usize, r: f64) -> Vec<f64> {
    let dx = L / n as f64;
    let m = n + 1;
    let tout: Vec<f64> = (0..m).map(|i| {
        let x = i as f64 * dx;
        if x < L / 2.0 { 20.0 } else if x > L / 2.0 { 27.0 } else { 23.5 }
    }).collect();
    let s = D * tau_of(r) / (dx * dx);                     // (lambda / dx)^2
    let mut a = vec![-s; m]; a[0] = 0.0; a[m - 1] = -2.0 * s;
    let mut c = vec![-s; m]; c[0] = -2.0 * s; c[m - 1] = 0.0;
    thomas(&a, &vec![1.0 + 2.0 * s; m], &c, &tout)
}

/// x where the increasing profile t crosses level v (linear interpolation).
fn cross(t: &[f64], x: &[f64], v: f64) -> f64 {
    let k = t.iter().position(|&ti| ti >= v).unwrap();
    x[k - 1] + (v - t[k - 1]) / (t[k] - t[k - 1]) * (x[k] - x[k - 1])
}
fn width(t: &[f64], x: &[f64]) -> f64 { cross(t, x, 26.3) - cross(t, x, 20.7) }

fn main() {
    let n = 200;
    let x: Vec<f64> = (0..=n).map(|i| i as f64 * L / n as f64).collect();
    let t = embryo_steady(n, 100.0);
    let lam = (D * tau_of(100.0)).sqrt();
    println!("tau {:.5} ms  lambda {:.4} um  width {:.4} um  formula {:.4} um",
             1000.0 * tau_of(100.0), lam, width(&t, &x), 2.0 * lam * 5f64.ln());
    let err = x.iter().zip(&t).map(|(xi, ti)| {
        let s = if *xi > L / 2.0 { 1.0 } else if *xi < L / 2.0 { -1.0 } else { 0.0 };
        (ti - (23.5 + 3.5 * s * (1.0 - (-(xi - L / 2.0).abs() / lam).exp()))).abs()
    }).fold(0.0, f64::max);
    println!("{:.6e}", err);                                  // formula ignores the tips
    let t200: Vec<f64> = [25, 50, 100, 200, 400].iter()
        .map(|&m| embryo_steady(m, 100.0)[m * 200 / 500]).collect();
    let d: Vec<f64> = t200.windows(2).map(|w| w[1] - w[0]).collect();
    println!("{:?}", d.windows(2).map(|w| format!("{:.2}", w[0] / w[1])).collect::<Vec<_>>());
    println!("{:?}", [50.0, 100.0, 150.0].iter()
        .map(|&r| format!("{:.2}", width(&embryo_steady(n, r), &x))).collect::<Vec<_>>());
}
