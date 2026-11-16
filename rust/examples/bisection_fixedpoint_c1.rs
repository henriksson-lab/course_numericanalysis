// C1: bisection, its step count, the floating-point floor, and a pole.
struct Bracket {
    root: f64,
    iter: usize,
    a: f64,
    b: f64,
}

// Result<..., String>: either Ok(bracket) or Err(message), like stop() in R
fn bisect(f: impl Fn(f64) -> f64, mut a: f64, mut b: f64, tol: f64, maxit: usize)
    -> Result<Bracket, String> {
    let mut fa = f(a);
    if fa.signum() == f(b).signum() {
        return Err("f(a) and f(b) must have opposite signs".to_string());
    }
    for k in 1..=maxit {
        let m = (a + b) / 2.0;
        let fm = f(m);
        if fm.signum() == fa.signum() { a = m; fa = fm; } else { b = m; }
        if b - a < tol * m.abs().max(1.0) {
            return Ok(Bracket { root: (a + b) / 2.0, iter: k, a, b });
        }
    }
    eprintln!("warning: maxit reached");
    Ok(Bracket { root: (a + b) / 2.0, iter: maxit, a, b })
}

fn main() {
    let sq = |x: f64| x * x - 3.0;
    for tol in [1e-6, 1e-10, 1e-14] {
        let out = bisect(sq, 1.0, 2.0, tol, 100).unwrap();
        let pred = (1.0 / (tol * 3f64.sqrt())).log2().ceil();
        println!("tol {tol:.0e}: {} steps (predicted {pred}), error {:.2e}",
                 out.iter, out.root - 3f64.sqrt());
    }
    let out = bisect(sq, 1.0, 2.0, 1e-17, 100).unwrap(); // warning: maxit reached
    println!("{} {:e}", out.iter, out.b - out.a);
    let pole = |x: f64| 1.0 / (x - 1.0);
    let out = bisect(pole, 0.0, 3.0, 1e-8, 100).unwrap();
    println!("{} {}", out.root, pole(out.root));
}
