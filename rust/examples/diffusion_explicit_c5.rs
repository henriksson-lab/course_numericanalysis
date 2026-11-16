// C5: no-flux walls: a dead-end chamber, and mass conservation with both ends closed.
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

// wall (no flux) at x = L via the ghost point c[N+1] = c[N-1]; optionally also at x = 0
fn ftcs_closed(c: &[f64], r: f64, nsteps: usize, left_closed: bool) -> Vec<f64> {
    let n = c.len();
    let mut c = c.to_vec();
    let mut old = c.clone();
    for _ in 0..nsteps {
        old.copy_from_slice(&c);
        for i in 1..n - 1 {
            c[i] = old[i] + r * (old[i - 1] - 2.0 * old[i] + old[i + 1]);
        }
        c[n - 1] = old[n - 1] + 2.0 * r * (old[n - 2] - old[n - 1]);
        if left_closed {
            c[0] = old[0] + 2.0 * r * (old[1] - old[0]);
        }
    }
    c
}

fn main() {
    let n = 50;
    let dx = L / n as f64;
    let x: Vec<f64> = (0..=n).map(|i| i as f64 * dx).collect();
    let (r, dt) = (0.4, 0.4 * dx * dx / D);
    let mut c0 = vec![0.0; n + 1];
    c0[0] = 1.0;
    // (a) open vs dead-end chamber at t = 64.4 s (gradient time of the open chamber)
    let nst = (64.4 / dt).round() as usize;
    let open = ftcs(&c0, r, nst);
    let dead = ftcs_closed(&c0, r, nst, false);
    for j in (0..=n).step_by(10) {
        println!("x = {:3.0} um   open {:.3}   dead end {:.3}", x[j], open[j], dead[j]);
    }
    // (b) time until the closed end reaches 0.95
    let mut c = c0.clone();
    let mut steps = 0;
    while c[n] < 0.95 {
        c = ftcs_closed(&c, r, 1, false);
        steps += 1;
    }
    let t = steps as f64 * dt;
    println!("closed end reaches 0.95 after {:.1} s; t D / L^2 = {:.3}", t, t * D / (L * L));
    // (c) both ends closed: a bolus of chemokine in the middle
    let b: Vec<f64> = x.iter().map(|&xi| if (xi - 250.0).abs() <= 20.0 { 1.0 } else { 0.0 }).collect();
    let w: Vec<f64> = (0..=n).map(|i| if i == 0 || i == n { 0.5 * dx } else { dx }).collect(); // trapezoid
    for t in [0.0, 5.0, 60.0] {
        let bt = ftcs_closed(&b, r, (t / dt).round() as usize, true);
        let trap: f64 = w.iter().zip(&bt).map(|(a, b)| a * b).sum();
        let plain: f64 = bt.iter().sum::<f64>() * dx;
        let mx = bt.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        println!("t = {:3.0} s   trapezoid amount {:.10}   plain sum*dx {:.6}   max {:.4}", t, trap, plain, mx);
    }
}
