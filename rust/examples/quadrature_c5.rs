// C5: what the quadrature crate (double-exponential rule) sees. Its points
// are packed towards both ends of the interval, unlike QUADPACK's.
use statrs::distribution::{Continuous, Normal};
use std::cell::RefCell;

fn main() {
    let dnorm = |x: f64| Normal::standard().pdf(x);
    // (a) record the points; RefCell lets the closure push to a vector
    let seen = RefCell::new(Vec::new());
    let spy = |x: f64| {
        seen.borrow_mut().push(x);
        dnorm(x)
    };
    let out = quadrature::integrate(spy, 0.0, 20000.0, 1e-8);
    let mut s = seen.borrow().clone();
    s.sort_by(|a, b| a.partial_cmp(b).unwrap());
    println!("{} with error estimate {:.1e}", out.integral, out.error_estimate);
    println!("{} points, smallest: {:.2e} {:.2e} {:.2e}", s.len(), s[0], s[1], s[2]);
    // (b) vary the upper limit (no infinite limits in this crate)
    for b in [20.0, 200.0, 2000.0, 5000.0, 20000.0] {
        let o = quadrature::integrate(dnorm, 0.0, b, 1e-8);
        println!("b = {:5}: {:.6} (error estimate {:.1e})", b, o.integral, o.error_estimate);
    }
    // (c) a peak inside the interval: the points miss it; split at the feature
    let shifted = |x: f64| dnorm(x - 5000.0);
    println!("{}", quadrature::integrate(shifted, 0.0, 20000.0, 1e-8).integral);
    let parts = quadrature::integrate(shifted, 0.0, 5000.0, 1e-8).integral
        + quadrature::integrate(shifted, 5000.0, 20000.0, 1e-8).integral;
    println!("{}", parts);
    // (d) a one-hour chromatogram with one peak at 340 s
    let peak = |t: f64| 1000.0 * (-(t - 340.0).powi(2) / (2.0 * 9.0)).exp();
    let o = quadrature::integrate(peak, 0.0, 3600.0, 1e-6);
    println!("{:e} with error estimate {:.1e}", o.integral, o.error_estimate);
    let o = quadrature::integrate(peak, 310.0, 370.0, 1e-6); // +-10 sigma
    println!("{} with error estimate {:.1e}", o.integral, o.error_estimate);
}
