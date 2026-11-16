// C5: integrate the standard normal density from 0 to 20000, 20 and infinity.
// R's integrate and scipy's quad are replaced by quadrature::integrate, a
// double-exponential rule. It packs its points densely near both ends of the
// interval, so here it finds the peak at 0 that R and SciPy miss.
use statrs::distribution::{Continuous, Normal};

fn main() {
    let dnorm = |x: f64| Normal::standard().pdf(x); // a closure, like R's dnorm
    for b in [20000.0, 20.0] {
        let out = quadrature::integrate(dnorm, 0.0, b, 1e-8);
        println!("0 to {}: {} with error estimate {:.1e}", b, out.integral, out.error_estimate);
    }
    // The crate has no infinite limits; substitute x = t / (1 - t), which maps
    // t in [0, 1) to x in [0, inf), with dx = dt / (1 - t)^2
    let g = |t: f64| dnorm(t / (1.0 - t)) / (1.0 - t).powi(2);
    let out = quadrature::integrate(g, 0.0, 1.0, 1e-8);
    println!("0 to inf: {} with error estimate {:.1e}", out.integral, out.error_estimate);
    // Move the peak away from the ends (a normal density centred at 5000):
    // now every point misses it, and the answer is 0 with error estimate 0.
    let out = quadrature::integrate(|x| dnorm(x - 5000.0), 0.0, 20000.0, 1e-8);
    println!("peak at 5000, 0 to 20000: {} with error estimate {:.1e}", out.integral, out.error_estimate);
}
