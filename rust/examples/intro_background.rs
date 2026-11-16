// Errors of the Taylor polynomials P1, P2, P3 for e^0.1.
fn main() {
    let x: f64 = 0.1; // f64 = a 64-bit floating-point number, like R's numeric
    let p = [1.0 + x, 1.0 + x + x * x / 2.0, 1.0 + x + x * x / 2.0 + x.powi(3) / 6.0];
    // .map applies the closure |pk| ... to each element, like sapply in R
    let err = p.map(|pk| x.exp() - pk);
    // {:.6e} prints in scientific notation with 6 decimals (Rust writes e-3, not e-03)
    for e in err {
        print!("{:.6e} ", e);
    }
    println!();
}
