fn main() {
    for k in 3..=17 {
        let mu = 10f64.powi(-k);
        let naive = 1.0 - (-mu).exp();
        let stable = -(-mu).exp_m1();           // exp_m1 is Rust's expm1
        let relerr = (naive - stable).abs() / stable;
        println!("{:.0e}  {:.1e}", mu, relerr);
    }
}
