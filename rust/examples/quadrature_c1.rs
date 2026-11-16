// C1: composite left Riemann, midpoint, trapezoid and Simpson rules on the
// AUC 0-6 h of the drug model; errors and the order (slope on log-log axes).
fn model(t: f64) -> f64 {
    500.0 / 50.0 * 1.5 / 1.3 * ((-0.2 * t).exp() - (-1.5 * t).exp())
}
fn left(f: fn(f64) -> f64, a: f64, b: f64, n: usize) -> f64 {
    let h = (b - a) / n as f64;
    h * (0..n).map(|i| f(a + i as f64 * h)).sum::<f64>()
}
fn mid(f: fn(f64) -> f64, a: f64, b: f64, n: usize) -> f64 {
    let h = (b - a) / n as f64;
    h * (0..n).map(|i| f(a + (i as f64 + 0.5) * h)).sum::<f64>()
}
fn trap(f: fn(f64) -> f64, a: f64, b: f64, n: usize) -> f64 {
    let h = (b - a) / n as f64;
    let inner: f64 = (1..n).map(|i| f(a + i as f64 * h)).sum();
    h * (inner + (f(a) + f(b)) / 2.0)
}
fn simpson(f: fn(f64) -> f64, a: f64, b: f64, n: usize) -> f64 {
    let h = (b - a) / n as f64;
    let s: f64 = (0..=n)
        .map(|i| {
            let w = if i == 0 || i == n { 1.0 } else if i % 2 == 1 { 4.0 } else { 2.0 };
            w * f(a + i as f64 * h)
        })
        .sum();
    h / 3.0 * s
}
// slope of log(e) against log(n) by least squares
fn slope(n: &[f64], e: &[f64]) -> f64 {
    let x: Vec<f64> = n.iter().map(|v| v.ln()).collect();
    let y: Vec<f64> = e.iter().map(|v| v.ln()).collect();
    let (mx, my) = (x.iter().sum::<f64>() / x.len() as f64, y.iter().sum::<f64>() / y.len() as f64);
    let sxy: f64 = x.iter().zip(&y).map(|(a, b)| (a - mx) * (b - my)).sum();
    let sxx: f64 = x.iter().map(|a| (a - mx).powi(2)).sum();
    sxy / sxx
}

fn main() {
    let exact = 32.6243601574044; // AUC 0-6 h
    let rules: [fn(fn(f64) -> f64, f64, f64, usize) -> f64; 4] = [left, mid, trap, simpson];
    let ns: Vec<usize> = (1..=10).map(|k| 1usize << k).collect(); // 2, 4, ..., 1024
    let err = |r: usize, n: usize| (rules[r](model, 0.0, 6.0, n) - exact).abs();
    for &n in &[4, 8, 16, 1024] {
        let row: Vec<String> = (0..4).map(|r| format!("{:.3e}", err(r, n))).collect();
        println!("{} {}", n, row.join(" "));
    }
    let big: Vec<usize> = ns.iter().copied().filter(|&n| n >= 32).collect();
    let nf: Vec<f64> = big.iter().map(|&n| n as f64).collect();
    let slopes: Vec<String> = (0..4)
        .map(|r| format!("{:.4}", slope(&nf, &big.iter().map(|&n| err(r, n)).collect::<Vec<_>>())))
        .collect();
    println!("slopes: {}", slopes.join(" "));
}
