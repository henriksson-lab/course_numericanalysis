// Background: a hand-written Euler method on the logistic growth curve.
const K: f64 = 1.2; // OD units
const N0: f64 = 0.05;
const R: f64 = 0.8; // per hour

fn f(_t: f64, n: f64) -> f64 {
    R * n * (1.0 - n / K) // logistic growth: f(t, y); `_t` = argument not used
}

fn exact(t: f64) -> f64 {
    K / (1.0 + (K / N0 - 1.0) * (-R * t).exp())
}

// returns (t, y) as two vectors; `impl Fn(f64, f64) -> f64` accepts any f(t, y)
fn euler(f: impl Fn(f64, f64) -> f64, t0: f64, y0: f64, h: f64, n: usize) -> (Vec<f64>, Vec<f64>) {
    let t: Vec<f64> = (0..=n).map(|i| t0 + i as f64 * h).collect();
    let mut y = vec![0.0; n + 1];
    y[0] = y0; // y[i] holds y_i
    for i in 0..n {
        y[i + 1] = y[i] + h * f(t[i], y[i]);
    }
    (t, y)
}

fn main() {
    let (t, y) = euler(f, 0.0, N0, 1.0, 5);
    for i in 0..t.len() {
        println!("{:3.0}  {:.5}  {:.5}", t[i], y[i], exact(t[i]));
    }
}
