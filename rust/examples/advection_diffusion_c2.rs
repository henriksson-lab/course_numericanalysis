// C2: a 1 s pulse carried to cells 2 mm downstream, explicit upwind scheme.
// Numerical diffusion (peak vs dx) and the stability condition C + 2r <= 1.
use statrs::function::erf::erfc;

const U: f64 = 200.0; // um/s
const D: f64 = 100.0; // um^2/s

struct Run { peak: f64, tpeak: f64, r: f64, maxabs: f64 }

fn pulse(dx: f64, cfl: f64, tend: f64) -> Run {
    let n = (3000.0 / dx).round() as usize + 1;
    let iobs = (2000.0 / dx).round() as usize;
    let dt = cfl * dx / U;
    let r = D * dt / (dx * dx);
    let mut c = vec![0.0; n];
    let mut old = c.clone();
    let nsteps = (tend / dt).round() as usize;
    let (mut peak, mut kpeak, mut big) = (f64::NEG_INFINITY, 0, 0.0f64);
    for k in 1..=nsteps {
        old.copy_from_slice(&c); // previous time level
        for i in 1..n - 1 {
            c[i] = old[i] - cfl * (old[i] - old[i - 1]) + r * (old[i - 1] - 2.0 * old[i] + old[i + 1]);
        }
        c[n - 1] = c[n - 2]; // outflow
        c[0] = if k as f64 * dt <= 1.0 { 1.0 } else { 0.0 }; // inlet: 1 s pulse
        if c[iobs] > peak { peak = c[iobs]; kpeak = k; }
        big = c.iter().fold(big, |m, v| m.max(v.abs()));
    }
    Run { peak, tpeak: kpeak as f64 * dt, r, maxabs: big }
}

fn step(x: f64, t: f64, d: f64) -> f64 {
    if t > 0.0 { 0.5 * erfc((x - U * t) / (2.0 * (d * t).sqrt())) } else { 0.0 }
}

// like C's %.3g for the numbers printed here
fn g3(x: f64) -> String {
    if x.abs() >= 1e4 {
        let s = format!("{:.2e}", x);
        let (m, e) = s.split_once('e').unwrap();
        let e: i32 = e.parse().unwrap();
        format!("{}e{}{:02}", m, if e < 0 { '-' } else { '+' }, e.abs())
    } else {
        let s = format!("{:.*}", (2 - x.abs().log10().floor() as i32).max(0) as usize, x);
        if s.contains('.') { s.trim_end_matches('0').trim_end_matches('.').to_string() } else { s }
    }
}

fn main() {
    for dd in [100.0, 600.0] {
        let peak = (0..=2200).map(|j| 5.0 + 0.005 * j as f64)
            .map(|t| step(2000.0, t, dd) - step(2000.0, t - 1.0, dd)).fold(0.0, f64::max);
        println!("long-channel formula, D = {:3}: peak {:.4}", dd, peak);
    }
    for dx in [10.0, 5.0, 2.0, 1.0] {
        let p = pulse(dx, 0.5, 16.0);
        println!("dx = {:4.1}  D_num = {:5.1}  peak {:.4} at t = {:.3} s", dx, U * dx * 0.5 / 2.0, p.peak, p.tpeak);
    }
    for cfl in [0.5, 0.9, 0.95, 1.0] {
        let p = pulse(10.0, cfl, 10.0);
        println!("C = {:.2}  r = {:.4}  C + 2r = {:.3}  max|c| = {}", cfl, p.r, cfl + 2.0 * p.r, g3(p.maxabs));
    }
}
