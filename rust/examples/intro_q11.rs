// Q11: the normal density integrated over limits that fit it.
// R's integrate is replaced by the quadrature crate (double-exponential rule).
use statrs::distribution::{Continuous, Normal};

fn main() {
    let out = quadrature::integrate(|x| Normal::standard().pdf(x), 0.0, 20.0, 1e-8);
    println!("{} with error estimate {:.1e}", out.integral, out.error_estimate);
}
