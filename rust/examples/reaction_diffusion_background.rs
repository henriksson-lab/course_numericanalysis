// Background: steady synthesis-diffusion-degradation gradient on a grid with
// ghost points, solved by Thomas, checked against the exact cosh formula;
// then one hour of IMEX Euler from c = 0.

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

/// IMEX Euler: diffusion implicit (one Thomas solve per step), reaction f explicit.
fn imex(c: &[f64], lo: &[f64], di: &[f64], up: &[f64], src: &[f64], dt: f64, nsteps: usize,
        f: &dyn Fn(&[f64]) -> Vec<f64>) -> Vec<f64> {
    let a: Vec<f64> = lo.iter().map(|v| dt * v).collect();
    let b: Vec<f64> = di.iter().map(|v| 1.0 + dt * v).collect();
    let cc: Vec<f64> = up.iter().map(|v| dt * v).collect();
    let mut c = c.to_vec();
    for _ in 0..nsteps {
        let fc = f(&c);
        let rhs: Vec<f64> = (0..c.len()).map(|i| c[i] + dt * (fc[i] + src[i])).collect();
        c = thomas(&a, &b, &cc, &rhs);
    }
    c
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
    let di: Vec<f64> = vec![2.0 * r; n + 1];
    let di_k: Vec<f64> = di.iter().map(|v| v + k).collect(); // + decay
    let mut src = vec![0.0; n + 1];
    src[0] = 2.0 * j / dx; // flux j enters at x = 0
    let cs = thomas(&lo, &di_k, &up, &src); // steady state: -D c'' + k c = src
    let err = x.iter().zip(&cs)
        .map(|(xi, ci)| (ci - j * lam / dd * ((l - xi) / lam).cosh() / (l / lam).sinh()).abs())
        .fold(0.0, f64::max);
    println!("{}", err);
    let total: f64 = cs.iter().sum::<f64>() - (cs[0] + cs[n]) / 2.0;
    println!("{}", k * dx * total); // total decay (trapezoid) = influx j
    let c60 = imex(&vec![0.0; n + 1], &lo, &di, &up, &src, 1.0, 60,
                   &|c: &[f64]| c.iter().map(|v| -k * v).collect());
    println!("{:.4} {:.4} {:.4}", c60[0], c60[10], c60[20]); // x = 0, 100, 200 µm at t = 60 min
}
