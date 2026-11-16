// Background: FTCS for the chemotaxis chamber (source at x = 0, sink at x = L).
use std::f64::consts::PI;

const D: f64 = 1000.0; // µm^2/s: small molecule, 1e-9 m^2/s
const L: f64 = 500.0; // µm

// FTCS with fixed (Dirichlet) ends: c[0] and c[n-1] never change.
// `old` keeps the previous time level, so every update uses old values only.
fn ftcs(c: &[f64], r: f64, nsteps: usize) -> Vec<f64> {
    let n = c.len();
    let mut c = c.to_vec();
    let mut old = c.clone();
    for _ in 0..nsteps {
        old.copy_from_slice(&c);
        for i in 1..n - 1 {
            c[i] = old[i] + r * (old[i - 1] - 2.0 * old[i] + old[i + 1]);
        }
    }
    c
}

// exact solution (Fourier series), chamber empty at t = 0
fn exact(x: f64, t: f64, l: f64, d: f64, nterm: usize) -> f64 {
    let mut s = 1.0 - x / l;
    for n in 1..=nterm {
        let k = n as f64 * PI / l;
        s -= 2.0 / (n as f64 * PI) * (k * x).sin() * (-k * k * d * t).exp();
    }
    s
}

fn main() {
    let n = 50; // 50 intervals of 10 µm, 51 grid points
    let dx = L / n as f64;
    let x: Vec<f64> = (0..=n).map(|i| i as f64 * dx).collect();
    let r = 0.4;
    let dt = r * dx * dx / D; // 0.04 s; FTCS needs r <= 1/2
    let mut c0 = vec![0.0; n + 1];
    c0[0] = 1.0; // source 1, channel empty, sink 0
    let c60 = ftcs(&c0, r, (60.0 / dt).round() as usize); // 1500 steps to t = 60 s
    for j in (5..n).step_by(10) {
        // x = 50, 150, ..., 450 µm
        println!("x = {:3.0} um   FTCS {:.4}   exact {:.4}", x[j], c60[j], exact(x[j], 60.0, L, D, 200));
    }
}
