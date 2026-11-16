// One-pass vs two-pass variance, and the shift check.
// Note: .sum() adds left to right in plain doubles. R's sum() uses extended
// precision and numpy's uses pairwise summation, so *wrong* one-pass results
// can differ between R, Python and Rust; the two-pass results agree.

fn mean(x: &[f64]) -> f64 {
    x.iter().sum::<f64>() / x.len() as f64
}

fn one_pass(x: &[f64]) -> f64 {
    let sq: Vec<f64> = x.iter().map(|v| v * v).collect();
    mean(&sq) - mean(x).powi(2)
}

fn two_pass(x: &[f64]) -> f64 {
    let m = mean(x);
    let dev2: Vec<f64> = x.iter().map(|v| (v - m).powi(2)).collect();
    mean(&dev2)
}

fn main() {
    // (a)
    for m in [1e4, 1e6, 1e8, 1e9] {
        // values m-3, ..., m+3, 100 of each; `i % 7` is the remainder
        let x: Vec<f64> = (0..700).map(|i| m + (i % 7) as f64 - 3.0).collect();
        println!("m = {:.0e}   one-pass = {}   two-pass = {}", m, one_pass(&x), two_pass(&x));
    }

    // (b)
    let x: Vec<f64> = [0.1, 0.2, 0.3].iter().map(|v| 7e8 + v).collect();
    println!("{}", one_pass(&x)); // one-pass
    println!("{}", two_pass(&x)); // two-pass

    // (c) shift check (Rust counts from 0, so x[0] is R's x[1])
    let x0: Vec<f64> = x.iter().map(|v| v - x[0]).collect();
    println!("{}", one_pass(&x0));
}
