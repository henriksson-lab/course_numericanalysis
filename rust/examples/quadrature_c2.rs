// C2: order drops for sqrt(x); exp(-x)/sqrt(x) by the midpoint rule, by
// substitution x = u^2, and by the quadrature crate.
fn trap(f: impl Fn(f64) -> f64, a: f64, b: f64, n: usize) -> f64 {
    let h = (b - a) / n as f64;
    let inner: f64 = (1..n).map(|i| f(a + i as f64 * h)).sum();
    h * (inner + (f(a) + f(b)) / 2.0)
}
fn simpson(f: impl Fn(f64) -> f64, a: f64, b: f64, n: usize) -> f64 {
    let h = (b - a) / n as f64;
    let s: f64 = (0..=n)
        .map(|i| {
            let w = if i == 0 || i == n { 1.0 } else if i % 2 == 1 { 4.0 } else { 2.0 };
            w * f(a + i as f64 * h)
        })
        .sum();
    h / 3.0 * s
}
fn mid(f: impl Fn(f64) -> f64, a: f64, b: f64, n: usize) -> f64 {
    let h = (b - a) / n as f64;
    h * (0..n).map(|i| f(a + (i as f64 + 0.5) * h)).sum::<f64>()
}
fn slope(n: &[f64], e: &[f64]) -> f64 {
    let x: Vec<f64> = n.iter().map(|v| v.ln()).collect();
    let y: Vec<f64> = e.iter().map(|v| v.ln()).collect();
    let (mx, my) = (x.iter().sum::<f64>() / x.len() as f64, y.iter().sum::<f64>() / y.len() as f64);
    let sxy: f64 = x.iter().zip(&y).map(|(a, b)| (a - mx) * (b - my)).sum();
    sxy / x.iter().map(|a| (a - mx).powi(2)).sum::<f64>()
}

fn main() {
    // (a) sqrt(x) on [0, 1], exact 2/3
    let ns: Vec<usize> = (4..=10).map(|k| 1usize << k).collect();
    let nf: Vec<f64> = ns.iter().map(|&n| n as f64).collect();
    let e_t: Vec<f64> = ns.iter().map(|&n| (trap(f64::sqrt, 0.0, 1.0, n) - 2.0 / 3.0).abs()).collect();
    let e_s: Vec<f64> = ns.iter().map(|&n| (simpson(f64::sqrt, 0.0, 1.0, n) - 2.0 / 3.0).abs()).collect();
    println!("{} {}", slope(&nf, &e_t), slope(&nf, &e_s));
    // (b) exp(-x) / sqrt(x) on (0, 1]; exact = sqrt(pi) erf(1)
    let f = |x: f64| (-x).exp() / x.sqrt();
    let exact = std::f64::consts::PI.sqrt() * statrs::function::erf::erf(1.0);
    println!("{}", exact);
    println!("{}", trap(f, 0.0, 1.0, 16)); // f(0) = inf
    for n in [16, 256, 4096] {
        println!("midpoint n = {}: {:e}", n, mid(f, 0.0, 1.0, n) - exact);
    }
    // (c) substitute x = u^2: integrand 2 exp(-u^2), smooth
    let g = |u: f64| 2.0 * (-u * u).exp();
    for n in [2, 4, 8, 16] {
        println!("Simpson n = {}: {:e}", n, simpson(g, 0.0, 1.0, n) - exact);
    }
    // the crate replaces the inf at x = 0 by 0 and packs points near the ends
    let out = quadrature::integrate(f, 0.0, 1.0, 1e-10);
    println!("{} with error estimate {:.1e}", out.integral, out.error_estimate);
}
