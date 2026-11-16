// Cancellation curve: (1 + d) - 1 for d = 1e-1, ..., 1e-15.
// There is no plot in Rust (no plotting crate): make the log-log plot of
// relative error against d, with the line eps / d, in R or Python.
fn main() {
    for k in 1..=15 {
        let d = 10f64.powi(-k);
        let computed = (1.0 + d) - 1.0;
        let relerr = (computed - d).abs() / d;
        // {:>5.0e}: scientific notation, 0 decimals, right-aligned in 5 characters
        println!("{:>5.0e}  {:>7.1e}  {:4.1}", d, relerr, -relerr.log10());
    }
}
