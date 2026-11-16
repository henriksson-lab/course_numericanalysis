// C1: Euler on the logistic curve: check B3, errors for four step sizes, local error.
const K: f64 = 1.2;
const N0: f64 = 0.05;
const R: f64 = 0.8;

fn f(_t: f64, n: f64) -> f64 {
    R * n * (1.0 - n / K)
}
fn exact(t: f64) -> f64 {
    K / (1.0 + (K / N0 - 1.0) * (-R * t).exp())
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
    // (a) check B3
    println!("{:.5?}", euler(f, 0.0, N0, 1.0, 3).1);
    println!("{:.5?}", euler(f, 0.0, N0, 0.5, 2).1);

    // (b) error at t = 5 and largest error on [0, 10]
    for h in [1.0, 0.5, 0.1, 0.01] {
        let n = (10.0 / h).round() as usize;
        let (t, y) = euler(f, 0.0, N0, h, n);
        let i5 = t.iter().position(|&ti| (ti - 5.0).abs() < 1e-9).unwrap(); // index of t = 5
        let maxerr = t.iter().zip(&y).map(|(&ti, &yi)| (yi - exact(ti)).abs()).fold(0.0, f64::max);
        println!("{:5.2}  {:.6}  {:10.3e}  {:10.3e}", h, y[i5], y[i5] - exact(5.0), maxerr);
    }
    // No plotting here: write t and y to a file and plot in R or Python.

    // (c) local error of one step from t = 0
    let d2 = R * f(0.0, N0) * (1.0 - 2.0 * N0 / K); // y''(0)
    for h in [1.0, 0.5, 0.25, 0.1, 0.01] {
        println!("{:5.2}  {:10.3e}  {:10.3e}", h, N0 + h * f(0.0, N0) - exact(h), h * h / 2.0 * d2);
    }
}
