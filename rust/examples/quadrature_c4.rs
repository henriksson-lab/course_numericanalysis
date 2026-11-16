// C4: adaptive Simpson, written by hand; evaluation counts, and a narrow
// chromatography peak that the first five points miss.
use std::cell::Cell;

fn model(t: f64) -> f64 {
    500.0 / 50.0 * 1.5 / 1.3 * ((-0.2 * t).exp() - (-1.5 * t).exp())
}
fn adapt(f: &dyn Fn(f64) -> f64, a: f64, b: f64, fa: f64, fm: f64, fb: f64, tol: f64) -> f64 {
    let m = (a + b) / 2.0;
    let (fl, fr) = (f((a + m) / 2.0), f((m + b) / 2.0));
    let s1 = (b - a) / 6.0 * (fa + 4.0 * fm + fb);
    let s2 = (b - a) / 12.0 * (fa + 4.0 * fl + 2.0 * fm + 4.0 * fr + fb);
    if (s2 - s1).abs() < 15.0 * tol {
        return s2 + (s2 - s1) / 15.0;
    }
    adapt(f, a, m, fa, fl, fm, tol / 2.0) + adapt(f, m, b, fm, fr, fb, tol / 2.0)
}
fn integ(f: &dyn Fn(f64) -> f64, a: f64, b: f64, tol: f64) -> f64 {
    adapt(f, a, b, f(a), f((a + b) / 2.0), f(b), tol)
}

fn main() {
    let i24 = 49.5252069010181; // AUC 0-24 h of the model
    // (a) count evaluations; Cell lets the closure change a counter
    let calls = Cell::new(0);
    let counted = |t: f64| {
        calls.set(calls.get() + 1);
        model(t)
    };
    for tol in [1e-4, 1e-6, 1e-8] {
        calls.set(0);
        let err = integ(&counted, 0.0, 24.0, tol) - i24;
        println!("tol {:.0e}: error {:9.2e}, {} evaluations", tol, err, calls.get());
    }
    // (b) a chromatogram: one peak at 340 s, sigma 3 s, height 1000, run 0-600 s
    let peak = |t: f64| 1000.0 * (-(t - 340.0).powi(2) / (2.0 * 9.0)).exp();
    println!("{}", 1000.0 * 3.0 * (2.0 * std::f64::consts::PI).sqrt()); // true area
    println!("{:e}", integ(&peak, 0.0, 600.0, 1e-6));
    // (c) fix: start from 60 pieces of 10 s
    let total: f64 = (0..60).map(|i| integ(&peak, 10.0 * i as f64, 10.0 * (i + 1) as f64, 1e-6 / 60.0)).sum();
    println!("{}", total);
}
