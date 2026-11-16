// C3: Romberg integration (repeated Richardson extrapolation of the trapezoid).
fn model(t: f64) -> f64 {
    500.0 / 50.0 * 1.5 / 1.3 * ((-0.2 * t).exp() - (-1.5 * t).exp())
}
fn trap(f: fn(f64) -> f64, a: f64, b: f64, n: usize) -> f64 {
    let h = (b - a) / n as f64;
    let inner: f64 = (1..n).map(|i| f(a + i as f64 * h)).sum();
    h * (inner + (f(a) + f(b)) / 2.0)
}
fn romberg(f: fn(f64) -> f64, a: f64, b: f64, k: usize) -> Vec<Vec<f64>> {
    let mut r = vec![vec![f64::NAN; k]; k];
    for i in 0..k {
        r[i][0] = trap(f, a, b, 1 << i); // 2^i intervals
        for j in 1..=i {
            r[i][j] = r[i][j - 1] + (r[i][j - 1] - r[i - 1][j - 1]) / (4f64.powi(j as i32) - 1.0);
        }
    }
    r
}

fn main() {
    let exact = 32.6243601574044;
    // (a) diagonal of the table: best estimate with 2^(k-1) + 1 points
    for k in 4..=9 {
        let r = romberg(model, 0.0, 6.0, k);
        println!("k = {}, {:4} points: error {:.2e}", k, (1 << (k - 1)) + 1, (r[k - 1][k - 1] - exact).abs());
    }
    // (b) sqrt(x) on [0, 1]
    let r = romberg(f64::sqrt, 0.0, 1.0, 8);
    let d: Vec<String> = (0..8).map(|i| format!("{:.1e}", (r[i][i] - 2.0 / 3.0).abs())).collect();
    println!("{}", d.join(" "));
}
