// C4: FTCS above and below r = 1/2, and the eigenvalues of the second-difference matrix.
use nalgebra::DMatrix;
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

// scientific notation as R and Python print it: 4.09e+15 rather than Rust's 4.09e15
fn sci(x: f64, digits: usize) -> String {
    let s = format!("{:.*e}", digits, x);
    let (m, e) = s.split_once('e').unwrap();
    let e: i32 = e.parse().unwrap();
    format!("{}e{}{:02}", m, if e < 0 { '-' } else { '+' }, e.abs())
}

fn main() {
    let n_int = 50; // intervals
    let dx = L / n_int as f64;
    let mut c0 = vec![0.0; n_int + 1];
    c0[0] = 1.0;
    // (a) run to t = 30 s with four values of r
    for rr in [0.45, 0.5, 0.52, 0.6] {
        let dt = rr * dx * dx / D;
        let nst = (30.0 / dt).round() as usize;
        let cr = ftcs(&c0, rr, nst);
        let mn = cr.iter().cloned().fold(f64::INFINITY, f64::min);
        let mx = cr.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        println!("r = {:.2}  dt = {:.4} s  steps = {:4}  min c = {:>10}  max c = {:>10}",
                 rr, dt, nst, sci(mn, 3), sci(mx, 3));
    }
    // (b) eigenvalues of the second-difference matrix (49 interior points)
    let n = n_int - 1;
    let a = DMatrix::from_fn(n, n, |i, j| {
        if i == j { -2.0 } else if i.abs_diff(j) == 1 { 1.0 } else { 0.0 }
    });
    let lam = a.symmetric_eigenvalues();
    let (lmin, lmax): (f64, f64) = (lam.min(), lam.max());
    let formula: Vec<f64> = (1..=n).map(|k| -4.0 * (k as f64 * PI / (2.0 * n_int as f64)).sin().powi(2)).collect();
    let fmin = formula.iter().cloned().fold(f64::INFINITY, f64::min);
    let fmax = formula.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    println!("eigen(): {:.6} to {:.6}   formula: {:.6} to {:.6}", lmin, lmax, fmin, fmax);
    let rates = [D / (dx * dx) * lmin, D / (dx * dx) * lmax]; // per second
    println!("rates {:.4} and {:.4} per s; Euler needs dt < {:.5} s (r < {:.4})",
             rates[0], rates[1], 2.0 / rates[0].abs(), 2.0 / rates[0].abs() * D / (dx * dx));
    println!("slowest time constant {:.2} s = L^2/(pi^2 D) = {:.2} s; stiffness ratio {:.0}",
             1.0 / rates[1].abs(), L * L / (PI * PI * D), rates[0] / rates[1]);
}
