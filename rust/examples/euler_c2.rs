// C2: order of Euler on the logistic curve, step halving, steps for an error of 1e-6.
const K: f64 = 1.2;
const N0: f64 = 0.05;
const R: f64 = 0.8;

fn f(_t: f64, n: f64) -> f64 {
    R * n * (1.0 - n / K)
}
fn exact(t: f64) -> f64 {
    K / (1.0 + (K / N0 - 1.0) * (-R * t).exp())
}
// only the final value is needed here, so no vectors are stored
fn euler_last(f: impl Fn(f64, f64) -> f64, t0: f64, y0: f64, h: f64, n: usize) -> f64 {
    let mut y = y0;
    for i in 0..n {
        y += h * f(t0 + i as f64 * h, y);
    }
    y
}

fn main() {
    // (a) error at t = 5 for h = 1, 1/2, ..., 1/1024
    let h: Vec<f64> = (0..=10).map(|k| 2f64.powi(-k)).collect();
    let y5: Vec<f64> = h.iter().map(|&hi| euler_last(f, 0.0, N0, hi, (5.0 / hi).round() as usize)).collect();
    let err: Vec<f64> = y5.iter().map(|y| y - exact(5.0)).collect();
    let errs: Vec<String> = err.iter().map(|e| format!("{:.3e}", e)).collect();
    println!("{}", errs.join(" "));
    let ratios: Vec<String> = (0..10).map(|i| format!("{:.3}", err[i] / err[i + 1])).collect();
    println!("{}", ratios.join(" ")); // ratio when h is halved
    // least-squares slope of log10|err| against log10 h
    let x: Vec<f64> = h.iter().map(|v| v.log10()).collect();
    let z: Vec<f64> = err.iter().map(|v| v.abs().log10()).collect();
    let (mx, mz) = (x.iter().sum::<f64>() / 11.0, z.iter().sum::<f64>() / 11.0);
    let sxz: f64 = x.iter().zip(&z).map(|(a, b)| (a - mx) * (b - mz)).sum();
    let sxx: f64 = x.iter().map(|a| (a - mx).powi(2)).sum();
    println!("slope {:.6}", sxz / sxx);

    // (b) step halving: y(h) - y(h/2) estimates the error of y(h/2)
    for i in [0, 3, 6, 9] {
        println!("{:.7}  {:.4e}  {:.4e}", h[i + 1], y5[i] - y5[i + 1], err[i + 1]);
    }

    // (c) steps needed for an error of 1e-6
    let c = err[10].abs() / h[10];
    let n = (5.0 / (1e-6 / c)).ceil() as usize;
    println!("C = {:.4}   n = {}", c, n);
    let yn = euler_last(f, 0.0, N0, 5.0 / n as f64, n);
    println!("error with n steps: {:.4e}", yn - exact(5.0));
}
