// Relative error and number of correct digits.
// f64 is Rust's double-precision number (same as R's numeric and numpy's float64).
fn rel_err(approx: f64, truth: f64) -> f64 {
    (approx - truth).abs() / truth.abs()
}

fn digits(approx: f64, truth: f64) -> f64 {
    -rel_err(approx, truth).log10()
}

fn main() {
    println!("{}", rel_err(51.0, 50.0));
    println!("{}", digits(3.1416, std::f64::consts::PI));
}
