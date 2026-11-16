// C4: error against iteration; slopes of log10(error) give the rates.
const R: f64 = 0.5;
const K: f64 = 100.0;
const H: f64 = 10.0;

fn g(n: f64) -> f64 { n * (R * (1.0 - n / K)).exp() - H }
fn g2(n: f64) -> f64 { K * (1.0 - (1.0 + H / n).ln() / R) }
fn f(n: f64) -> f64 { g(n) - n }
fn dg(n: f64) -> f64 { (R * (1.0 - n / K)).exp() * (1.0 - R * n / K) }

// least-squares slope of log10(e_k) against k, using only the straight part
fn slope(e: &[f64]) -> f64 {
    let pts: Vec<(f64, f64)> = e.iter().enumerate()
        .filter(|&(_, &v)| v < 1e-2 && v > 1e-12)
        .map(|(i, &v)| ((i + 1) as f64, v.log10()))
        .collect();
    let n = pts.len() as f64;
    let mx = pts.iter().map(|p| p.0).sum::<f64>() / n;
    let my = pts.iter().map(|p| p.1).sum::<f64>() / n;
    let sxy: f64 = pts.iter().map(|p| (p.0 - mx) * (p.1 - my)).sum();
    let sxx: f64 = pts.iter().map(|p| (p.0 - mx).powi(2)).sum();
    sxy / sxx
}

fn main() {
    // reference root by bisection to full precision
    let (mut a, mut b) = (50.0, 100.0);
    for _ in 0..100 {
        let m: f64 = (a + b) / 2.0;
        if f(m).signum() == f(a).signum() { a = m } else { b = m }
    }
    let hi = (a + b) / 2.0;

    let n = 60;
    let (mut a, mut b) = (50.0, 100.0);
    let mut e_bis = vec![0.0; n];
    let (mut e_g, mut e_g2) = (vec![0.0; n], vec![0.0; n]);
    let (mut x, mut y) = (30.0, 30.0);
    for k in 0..n {
        let m = (a + b) / 2.0;
        if f(m).signum() == f(a).signum() { a = m } else { b = m }
        e_bis[k] = (b - a) / 2.0; // error bound of the midpoint
        x = g(x);
        y = g2(y);
        e_g[k] = (x - hi).abs();
        e_g2[k] = (y - hi).abs();
    }
    println!("{:.4} {:.4} {:.4}", slope(&e_bis), slope(&e_g), slope(&e_g2));
    let dg2 = K * H / (R * hi * (hi + H)); // g2'(N*)
    println!("{:.4} {:.4} {:.4}", 0.5f64.log10(), dg(hi).abs().log10(), dg2.log10());
    // No plotting here: write e_bis, e_g, e_g2 to a file and plot in R or Python.
}
