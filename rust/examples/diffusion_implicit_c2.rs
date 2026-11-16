// C2: verification with an exact sine mode: observed order in time and space.

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

/// One backward Euler step; c includes the two (fixed) boundary values.
fn be_step(c: &mut [f64], r: f64) {
    let m = c.len() - 2;
    let mut d = c[1..=m].to_vec();
    d[0] += r * c[0];
    d[m - 1] += r * c[m + 1];
    let x = thomas(&vec![-r; m], &vec![1.0 + 2.0 * r; m], &vec![-r; m], &d);
    c[1..=m].copy_from_slice(&x);
}

/// One Crank-Nicolson step.
fn cn_step(c: &mut [f64], r: f64) {
    let m = c.len() - 2;
    let mut d: Vec<f64> = (1..=m).map(|i| c[i] + r / 2.0 * (c[i - 1] - 2.0 * c[i] + c[i + 1])).collect();
    d[0] += r / 2.0 * c[0];
    d[m - 1] += r / 2.0 * c[m + 1];
    let x = thomas(&vec![-r / 2.0; m], &vec![1.0 + r; m], &vec![-r / 2.0; m], &d);
    c[1..=m].copy_from_slice(&x);
}

const DD: f64 = 1000.0;
const L: f64 = 500.0;

fn sine_err(n: usize, dt: f64, step: fn(&mut [f64], f64), tend: f64) -> f64 {
    let pi = std::f64::consts::PI;
    let dx = L / n as f64;
    let x: Vec<f64> = (0..=n).map(|i| i as f64 * dx).collect();
    let mut c: Vec<f64> = x.iter().map(|xi| (pi * xi / L).sin()).collect();
    for _ in 0..(tend / dt).round() as usize { step(&mut c, DD * dt / (dx * dx)); }
    let amp = (-DD * (pi / L).powi(2) * tend).exp();
    x.iter().zip(&c).map(|(xi, ci)| (ci - amp * (pi * xi / L).sin()).abs()).fold(0.0, f64::max)
}

fn order(e: &[f64]) -> Vec<String> {
    e.windows(2).map(|w| format!("{:.2}", (w[0] / w[1]).log2())).collect()
}

fn main() {
    let dts: Vec<f64> = (0..7).map(|k| 6.0 / 2f64.powi(k)).collect();
    let e_be: Vec<f64> = dts.iter().map(|&h| sine_err(400, h, be_step, 30.0)).collect();
    let e_cn: Vec<f64> = dts.iter().map(|&h| sine_err(400, h, cn_step, 30.0)).collect();
    let e_cn50: Vec<f64> = dts.iter().map(|&h| sine_err(50, h, cn_step, 30.0)).collect();
    for k in 0..7 {
        println!("{:9.4}  {:9.3e}  {:9.3e}  {:9.3e}", dts[k], e_be[k], e_cn[k], e_cn50[k]);
    }
    println!("{:?}", order(&e_be));
    println!("{:?}", order(&e_cn));
    let e_dx: Vec<f64> = [10, 20, 40, 80, 160].iter().map(|&n| sine_err(n, 0.01, cn_step, 30.0)).collect();
    println!("{:?}", order(&e_dx));
}
