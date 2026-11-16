// C1: FTCS against the exact solution; the straight line is a fixed point.
use std::f64::consts::PI;

const D: f64 = 1000.0; // µm^2/s
const L: f64 = 500.0; // µm

fn ftcs(c: &[f64], r: f64, nsteps: usize) -> Vec<f64> {
    let n = c.len();
    let mut c = c.to_vec();
    let mut old = c.clone();
    for _ in 0..nsteps {
        old.copy_from_slice(&c); // previous time level
        for i in 1..n - 1 {
            c[i] = old[i] + r * (old[i - 1] - 2.0 * old[i] + old[i + 1]);
        }
    }
    c
}

fn exact(x: f64, t: f64, l: f64, d: f64, nterm: usize) -> f64 {
    let mut s = 1.0 - x / l;
    for n in 1..=nterm {
        let k = n as f64 * PI / l;
        s -= 2.0 / (n as f64 * PI) * (k * x).sin() * (-k * k * d * t).exp();
    }
    s
}

// scientific notation as R and Python print it: 5.11e-04 rather than Rust's 5.11e-4
fn sci(x: f64, digits: usize) -> String {
    let s = format!("{:.*e}", digits, x);
    let (m, e) = s.split_once('e').unwrap();
    let e: i32 = e.parse().unwrap();
    format!("{}e{}{:02}", m, if e < 0 { '-' } else { '+' }, e.abs())
}

fn maxabs(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(u, v)| (u - v).abs()).fold(0.0, f64::max)
}

fn main() {
    let n = 50;
    let dx = L / n as f64;
    let x: Vec<f64> = (0..=n).map(|i| i as f64 * dx).collect();
    let (r, dt) = (0.4, 0.4 * dx * dx / D);
    let mut c0 = vec![0.0; n + 1];
    c0[0] = 1.0;
    let lin: Vec<f64> = x.iter().map(|&xi| 1.0 - xi / L).collect();
    // (a) FTCS against the exact solution, and distance from the straight line
    for t in [10.0, 30.0, 60.0, 120.0, 300.0] {
        let ct = ftcs(&c0, r, (t / dt).round() as usize);
        let ex: Vec<f64> = x.iter().map(|&xi| exact(xi, t, L, D, 200)).collect();
        println!("t = {:3.0} s   max|FTCS - exact| = {}   max|c - line| = {}",
                 t, sci(maxabs(&ct, &ex), 2), sci(maxabs(&ct, &lin), 2));
    }
    // (c) the straight line is a fixed point of FTCS
    println!("one step from the line moves it by {}", sci(maxabs(&ftcs(&lin, r, 1), &lin), 1));
}
