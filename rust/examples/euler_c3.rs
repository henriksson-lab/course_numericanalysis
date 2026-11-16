// C3: instability of Euler: drug clearance and logistic growth with large steps.
const K: f64 = 1.2;
const N0: f64 = 0.05;
const R: f64 = 0.8;

fn f(_t: f64, n: f64) -> f64 {
    R * n * (1.0 - n / K)
}
fn euler(f: impl Fn(f64, f64) -> f64, t0: f64, y0: f64, h: f64, n: usize) -> (Vec<f64>, Vec<f64>) {
    let t: Vec<f64> = (0..=n).map(|i| t0 + i as f64 * h).collect();
    let mut y = vec![0.0; n + 1];
    y[0] = y0;
    for i in 0..n {
        y[i + 1] = y[i] + h * f(t[i], y[i]);
    }
    (t, y)
}

fn main() {
    // (a) drug clearance C' = -kC, ten steps of each size
    let k = 0.8;
    let fc = |_t: f64, c: f64| -k * c; // a closure: a function that can use k
    for h in [0.5, 2.0, 2.4, 2.6, 3.0] {
        let c10 = euler(fc, 0.0, 10.0, h, 10).1[10];
        println!("h = {:.1}  factor {:5.2}  C_10 = {:10.3e}  exact {:9.3e}",
                 h, 1.0 - k * h, c10, 10.0 * (-k * 10.0 * h).exp());
    }
    // (b) logistic growth, 100 steps: the last four values
    for h in [2.0, 2.4, 2.6, 3.0, 1.5] {
        let y = euler(f, 0.0, N0, h, 100).1;
        let last: Vec<String> = y[97..].iter().map(|v| format!("{:.4}", v)).collect();
        println!("h = {:.1}  rh = {:.2} : {}", h, R * h, last.join(" "));
    }
}
