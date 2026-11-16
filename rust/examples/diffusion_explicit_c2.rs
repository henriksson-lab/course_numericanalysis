// C2: time until the linear gradient is established, for two molecules and two lengths.
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

fn t_grad(d: f64, l: f64, tol: f64, n: usize, r: f64) -> f64 {
    let dx = l / n as f64;
    let dt = r * dx * dx / d;
    let lin: Vec<f64> = (0..=n).map(|i| 1.0 - i as f64 * dx / l).collect();
    let mut c = vec![0.0; n + 1];
    c[0] = 1.0;
    let mut steps = 0;
    while c.iter().zip(&lin).map(|(u, v)| (u - v).abs()).fold(0.0, f64::max) > tol {
        c = ftcs(&c, r, 1);
        steps += 1;
    }
    steps as f64 * dt
}

fn main() {
    for tol in [0.05, 0.01] {
        for d in [1000.0, 100.0] {
            for l in [500.0, 1000.0] {
                let tg = t_grad(d, l, tol, 50, 0.4);
                println!("tol {:.2}  D = {:4.0}  L = {:4.0}   t = {:7.1} s = {:5.1} min   t D / L^2 = {:.4}",
                         tol, d, l, tg, tg / 60.0, tg * d / (l * l));
            }
        }
    }
}
