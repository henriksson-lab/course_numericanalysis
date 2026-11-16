// C4: Newton with a wrong derivative converges only linearly.
fn f(c: f64) -> f64 {
    0.3 * c / (0.1 + c) + 0.7 * c.powi(2) / (25.0 + c.powi(2)) - 0.5
}
fn df(c: f64) -> f64 {
    0.03 / (0.1 + c).powi(2) + 35.0 * c / (25.0 + c.powi(2)).powi(2)
}
fn df_bad(c: f64) -> f64 {
    0.03 / (0.1 + c).powi(2) + 50.0 * c / (25.0 + c.powi(2)).powi(2) // weight 0.7 forgotten
}

fn main() {
    let c50 = 3.261305265055212;
    let mut x: f64 = 1.0;
    let mut e = Vec::new();
    for _ in 0..30 {
        x -= f(x) / df_bad(x);
        e.push((x - c50).abs());
    }
    let shown: Vec<String> = [0, 1, 2, 3, 4, 5, 9, 19, 29].iter().map(|&i| format!("{:.2e}", e[i])).collect();
    println!("{:?}", shown);
    let ratios: Vec<String> = (0..11).map(|i| format!("{:.3}", e[i + 1] / e[i])).collect();
    println!("{:?}", ratios); // ratio e_{k+1}/e_k
    println!("{}", 1.0 - df(c50) / df_bad(c50)); // predicted ratio
}
