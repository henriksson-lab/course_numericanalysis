// C8: normalising a posterior density (zero-truncated Poisson UMI rate),
// with the quadrature crate. Infinite limits by x = t / (1 - t).
fn loglik(l: f64, m: f64) -> f64 {
    m * (173.0 * l.ln() - 100.0 * l - 100.0 * (-(-l).exp()).ln_1p())
}
// integrate f over [0, b], or [0, inf) when b is infinite
fn integ(f: impl Fn(f64) -> f64, b: f64, tol: f64) -> f64 {
    if b.is_infinite() {
        let g = |t: f64| f(t / (1.0 - t)) / (1.0 - t).powi(2);
        quadrature::integrate(g, 0.0, 1.0, tol).integral
    } else {
        quadrature::integrate(f, 0.0, b, tol).integral
    }
}

fn main() {
    let (lhat, se) = (1.2184695258301268, 0.1325); // from Newton and secant methods
    // (a) raw likelihood, target absolute error 1e-8 (like scipy's epsabs)
    let l = |x: f64| loglik(x, 1.0).exp();
    for b in [10.0, 100.0, f64::INFINITY] {
        println!("raw, 0 to {}: {:e}", b, integ(l, b, 1e-8));
    }
    // (b) scaled to a maximum of 1
    let g = |x: f64| (loglik(x, 1.0) - loglik(lhat, 1.0)).exp();
    let z = integ(g, f64::INFINITY, 1e-10);
    println!("{}", z);
    println!("{}", loglik(lhat, 1.0) + z.ln());
    println!("{}", integ(|x| x * g(x) / z, f64::INFINITY, 1e-10)); // posterior mean
    // (c) 10 000 cells
    let g2 = |x: f64| (loglik(x, 100.0) - loglik(lhat, 100.0)).exp();
    for b in [10.0, 100.0, f64::INFINITY] {
        println!("10 000 cells, 0 to {}: {:e}", b, integ(g2, b, 1e-10));
    }
    let se2 = se / 10.0;
    let o = quadrature::integrate(g2, lhat - 10.0 * se2, lhat + 10.0 * se2, 1e-10);
    println!("{}", o.integral);
}
