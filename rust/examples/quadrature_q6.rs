// Q6: a 3 s wide chromatography peak in a one-hour run, with the quadrature
// crate (double-exponential rule), over the whole run and over the peak.
fn main() {
    let peak = |t: f64| 1000.0 * (-(t - 340.0).powi(2) / (2.0 * 9.0)).exp();
    let o = quadrature::integrate(peak, 0.0, 3600.0, 1e-6);
    println!("{:e} with error estimate {:.1e}", o.integral, o.error_estimate);
    let o = quadrature::integrate(peak, 310.0, 370.0, 1e-6); // the peak +- 10 sigma
    println!("{} with error estimate {:.1e}", o.integral, o.error_estimate);
}
