// C5: Euler on the 24-hour oscillator x'' = -w^2 x, written as a first-order system.
use std::f64::consts::PI;

fn euler_sys(f: impl Fn(f64, &[f64]) -> Vec<f64>, t0: f64, y0: &[f64], h: f64, n: usize) -> Vec<Vec<f64>> {
    let mut ys = vec![y0.to_vec()];
    for i in 0..n {
        let y = &ys[i];
        let slope = f(t0 + i as f64 * h, y);
        let next: Vec<f64> = y.iter().zip(&slope).map(|(yj, fj)| yj + h * fj).collect();
        ys.push(next);
    }
    ys
}

const W: f64 = 2.0 * PI / 24.0; // one cycle per 24 hours

fn osc(_t: f64, y: &[f64]) -> Vec<f64> {
    vec![y[1], -W * W * y[0]] // y = (x, x')
}

// (a) amplitude after 10 days
fn amp(h: f64) -> f64 {
    let n = (240.0 / h).round() as usize;
    let ys = euler_sys(osc, 0.0, &[1.0, 0.0], h, n);
    (ys[n][0].powi(2) + (ys[n][1] / W).powi(2)).sqrt()
}

fn main() {
    for h in [1.0, 0.1, 0.01] {
        let predicted = (1.0 + W * W * h * h).powf(240.0 / h / 2.0);
        println!("h = {:4.2}  amplitude {:10.4}  predicted {:10.4}", h, amp(h), predicted);
    }
    // (b) phase portrait: write ys (h = 1, 48 steps) to a file and plot in R or Python.

    // (c) largest h with amplitude error below 1 % after 10 days
    let hmax = 1.01f64.ln() / (120.0 * W * W);
    let n = (240.0 / hmax).round();
    println!("h = {:.5} hours, {} steps; amplitude {:.5}", hmax, n, amp(240.0 / n));
}
