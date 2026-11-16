fn main() {
    // {:e} prints scientific notation; overflow gives inf, 0/0 gives NaN, no error or warning
    println!("{:e}", 1e300 * 1e10 / 1e20);
    println!("{:e}", 1e300 / 1e20 * 1e10);
    println!("{}", 1e-3f64.powi(2000) / 2e-3f64.powi(2000));
    println!("{}", 2000.0 * (1e-3f64.ln() - 2e-3f64.ln())); // the log of the ratio is fine
}
