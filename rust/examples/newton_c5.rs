// C5: a safeguarded Newton that never leaves its sign-change bracket.
fn f(c: f64) -> f64 {
    0.3 * c / (0.1 + c) + 0.7 * c.powi(2) / (25.0 + c.powi(2)) - 0.5
}
fn df(c: f64) -> f64 {
    0.03 / (0.1 + c).powi(2) + 35.0 * c / (25.0 + c.powi(2)).powi(2)
}
fn g(u: f64) -> f64 {
    1.0 / (1.0 + (-u).exp()) - 0.5 // Michaelis-Menten, log substrate
}
fn dg(u: f64) -> f64 {
    (-u).exp() / (1.0 + (-u).exp()).powi(2)
}

// returns (root, iterations, bisection steps)
fn newton_safe(f: impl Fn(f64) -> f64, df: impl Fn(f64) -> f64,
               mut a: f64, mut b: f64, x0: f64, tol: f64, maxit: usize) -> (f64, usize, usize) {
    let mut fa = f(a);
    assert!(fa * f(b) <= 0.0, "f(a) and f(b) must have opposite signs");
    let (mut x, mut old_step, mut nbis) = (x0, b - a, 0);
    for k in 1..=maxit {
        let fx = f(x);
        if fx == 0.0 {
            return (x, k, nbis);
        }
        if fx.signum() == fa.signum() { a = x; fa = fx; } else { b = x; } // keep the sign change
        let mut xnew = x - fx / df(x); // Newton proposal
        if !xnew.is_finite() || xnew <= a || xnew >= b || (xnew - x).abs() > old_step.abs() / 2.0 {
            xnew = (a + b) / 2.0; // outside or too slow: bisect
            nbis += 1;
        }
        old_step = xnew - x;
        if old_step.abs() <= tol * xnew.abs() {
            return (xnew, k, nbis);
        }
        x = xnew;
    }
    eprintln!("warning: no convergence");
    (x, maxit, nbis)
}

fn show(r: (f64, usize, usize)) {
    println!("root {:.6e}  iter {}  bisections {}", r.0, r.1, r.2);
}

fn main() {
    show(newton_safe(f, df, 0.01, 100.0, 10.0, 1e-12, 100));
    show(newton_safe(g, dg, -5.0, 5.0, 2.3, 1e-12, 100));
    show(newton_safe(g, dg, -5.0, 5.0, 2.1773189849653, 1e-12, 100));
}
