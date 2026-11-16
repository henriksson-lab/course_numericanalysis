// One-pass vs two-pass population variance.
// Note: .sum() adds left to right in plain doubles. R's sum() uses extended
// precision and numpy's uses pairwise summation, so a *wrong* one-pass result
// can differ between the three languages; the two-pass result agrees.
fn main() {
    let x: Vec<f64> = [1.0, 2.0, 3.0, 4.0].iter().map(|v| 1e8 + v).collect();
    let n = x.len() as f64; // `as f64` converts the integer length to a float
    let mean = x.iter().sum::<f64>() / n;

    let one = x.iter().map(|v| v * v).sum::<f64>() / n - mean * mean;
    println!("{}", one); // one-pass
    let two = x.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / n;
    println!("{}", two); // two-pass
}
