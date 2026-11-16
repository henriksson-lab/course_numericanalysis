// Negative variance from the one-pass formula vs a two-pass (centred) variance.
// Note: .sum() adds left to right in plain doubles (R: extended precision,
// numpy: pairwise), so a wrong one-pass value can differ between languages.
fn main() {
    let x: Vec<f64> = [0.1, 0.2, 0.3].iter().map(|v| 7e8 + v).collect();
    let n = x.len() as f64;
    let mean = x.iter().sum::<f64>() / n;
    let mean_sq = x.iter().map(|v| v * v).sum::<f64>() / n;
    println!("{}", mean_sq - mean * mean); // one-pass: impossible
    // two-pass, population scaling (divides by n), as in C3
    let var = x.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / n;
    println!("{}", var);
    // The statrs crate's population_variance() (a running Welford update)
    // gives 0.006666655937834776 here: correct to about 6 digits, far better
    // than one-pass but less accurate than two-pass.
}
