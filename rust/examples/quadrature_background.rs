// Background: trapezoid rule on the blood samples (hand-written) and the
// quadrature crate's integrate (double-exponential rule) on the model.
fn auc_trap(t: &[f64], y: &[f64]) -> f64 {
    (1..t.len()).map(|i| (t[i] - t[i - 1]) * (y[i - 1] + y[i]) / 2.0).sum()
}
fn model(t: f64) -> f64 {
    500.0 / 50.0 * 1.5 / 1.3 * ((-0.2 * t).exp() - (-1.5 * t).exp())
}

fn main() {
    let t_h = [0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 24.0]; // hours
    let conc = [0.0, 4.99, 6.87, 7.16, 6.20, 5.16, 3.47, 2.33, 1.05, 0.09]; // mg/L
    println!("{}", auc_trap(&t_h, &conc));
    let out = quadrature::integrate(model, 0.0, 24.0, 1e-10); // target absolute error
    println!("{}", out.integral);
    println!("{:e}", out.error_estimate);
}
