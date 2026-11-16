// C3: grid refinement with r = 0.4 fixed: error at t = 30 s against the exact solution.
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

// scientific notation as R and Python print it: 1.11e-03 rather than Rust's 1.11e-3
fn sci(x: f64, digits: usize) -> String {
    let s = format!("{:.*e}", digits, x);
    let (m, e) = s.split_once('e').unwrap();
    let e: i32 = e.parse().unwrap();
    format!("{}e{}{:02}", m, if e < 0 { '-' } else { '+' }, e.abs())
}

fn main() {
    let ns = [20usize, 40, 80, 160, 320];
    let mut err = Vec::new();
    for &n in &ns {
        let dx = L / n as f64;
        let dt = 0.4 * dx * dx / D;
        let nst = (30.0 / dt).round() as usize; // steps to t = 30 s
        let mut c = vec![0.0; n + 1];
        c[0] = 1.0;
        let c = ftcs(&c, 0.4, nst);
        let e = (0..=n)
            .map(|i| (c[i] - exact(i as f64 * dx, 30.0, L, D, 2000)).abs())
            .fold(0.0, f64::max);
        err.push(e);
        println!("N = {:3}  dx = {:6.3} um  dt = {:8.5} s  steps = {:5}  work = {:8}  error = {}",
                 n, dx, dt, nst, nst * (n - 1), sci(e, 3));
    }
    let ratios: Vec<String> = err.windows(2).map(|w| format!("{:.2}", w[0] / w[1])).collect();
    println!("ratios: {}", ratios.join(" "));
    // least-squares slope of log10(error) against log10(dx)
    let lx: Vec<f64> = ns.iter().map(|&n| (L / n as f64).log10()).collect();
    let ly: Vec<f64> = err.iter().map(|e| e.log10()).collect();
    let (mx, my) = (lx.iter().sum::<f64>() / 5.0, ly.iter().sum::<f64>() / 5.0);
    let sxy: f64 = lx.iter().zip(&ly).map(|(a, b)| (a - mx) * (b - my)).sum();
    let sxx: f64 = lx.iter().map(|a| (a - mx) * (a - mx)).sum();
    println!("slope: {:.3}", sxy / sxx);
}
