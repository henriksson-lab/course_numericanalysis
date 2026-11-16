// C3: fixed-point iteration with three different g, rate and error estimates.
const R: f64 = 0.5;
const K: f64 = 100.0;
const H: f64 = 10.0;

fn g(n: f64) -> f64 { n * (R * (1.0 - n / K)).exp() - H }
fn g2(n: f64) -> f64 { K * (1.0 - (1.0 + H / n).ln() / R) }
fn g3(n: f64) -> f64 { H / ((R * (1.0 - n / K)).exp() - 1.0) }
fn f(n: f64) -> f64 { g(n) - n }

fn bisect(f: impl Fn(f64) -> f64, mut a: f64, mut b: f64, tol: f64) -> f64 {
    let mut fa = f(a);
    for _ in 0..200 {
        let m = (a + b) / 2.0;
        let fm = f(m);
        if fm.signum() == fa.signum() { a = m; fa = fm; } else { b = m; }
        if b - a < tol * m.abs().max(1.0) { break; }
    }
    (a + b) / 2.0
}

// returns (converged?, all iterates x0, x1, ...)
fn fixed_point(g: fn(f64) -> f64, x0: f64, tol: f64, maxit: usize) -> (bool, Vec<f64>) {
    let mut x = vec![x0];
    for _ in 0..maxit {
        let next = g(*x.last().unwrap());
        x.push(next);
        if !next.is_finite() { break; }
        let n = x.len();
        if (x[n - 1] - x[n - 2]).abs() < tol { return (true, x); }
    }
    eprintln!("warning: no convergence");
    (false, x)
}

fn main() {
    let hi = bisect(f, 50.0, 100.0, 1e-15);
    let lo = bisect(f, 0.0, 50.0, 1e-15);
    let runs: [(&str, fn(f64) -> f64, f64); 4] =
        [("g ", g, 30.0), ("g ", g, 21.0), ("g2", g2, 30.0), ("g3", g3, 30.0)];
    for (name, gf, x0) in runs {
        let (_, x) = fixed_point(gf, x0, 1e-10, 1000);
        let d: Vec<f64> = x.windows(2).map(|w| w[1] - w[0]).collect(); // steps
        let n = d.len();
        let l = d[n - 1] / d[n - 2];
        println!("{name} from {x0:2}: {:3} steps, root {:.10}, L = {l:.4}, est. error {:.1e}",
                 n, x[n], d[n - 1].abs() * l / (1.0 - l));
    }
    println!("{hi} {lo}");
    let (_, x) = fixed_point(g, 30.0, 1e-10, 1000);
    println!("{:e}", x.last().unwrap() - hi);
    let (_, x) = fixed_point(g, 20.0, 1e-10, 1000);       // warning: no convergence
    println!("{:?}", &x[x.len() - 4..]);
}
