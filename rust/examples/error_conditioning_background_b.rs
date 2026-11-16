// Condition number estimated by nudging the input by a relative amount h.
// `f: impl Fn(f64) -> f64` accepts any function or closure from f64 to f64.
fn kappa_num(f: impl Fn(f64) -> f64, x: f64, h: f64) -> f64 {
    let rel_out = (f(x * (1.0 + h)) - f(x)).abs() / f(x).abs();
    rel_out / h // relative change out / relative change in
}

fn main() {
    // Rust has no default arguments, so h = 1e-7 is passed explicitly.
    println!("{}", kappa_num(f64::ln, 1.001, 1e-7)); // f64::ln is the natural log
}
