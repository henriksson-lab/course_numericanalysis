// C3: the secant method; estimate the order of convergence from the errors.
const N: f64 = 100.0; // cells
const SX: f64 = 173.0; // total UMIs

fn f(c: f64) -> f64 {
    0.3 * c / (0.1 + c) + 0.7 * c.powi(2) / (25.0 + c.powi(2)) - 0.5
}
fn score(l: f64) -> f64 {
    SX / l - N / (-(-l).exp_m1())
}
fn dscore(l: f64) -> f64 {
    -SX / l.powi(2) + N * (-l).exp() / (-l).exp_m1().powi(2)
}

fn secant(f: impl Fn(f64) -> f64, mut x0: f64, mut x1: f64, tol: f64, maxit: usize) -> (f64, usize) {
    let (mut f0, mut f1) = (f(x0), f(x1));
    let mut k = 0;
    while k < maxit {
        k += 1;
        if f1 == f0 {
            break; // flat chord: stop
        }
        let step = f1 * (x1 - x0) / (f1 - f0);
        x0 = x1;
        f0 = f1;
        x1 -= step;
        f1 = f(x1);
        if step.abs() <= tol * x1.abs() {
            return (x1, k);
        }
    }
    eprintln!("warning: no convergence after {k} steps");
    (x1, k)
}

fn order_est(e: &[f64]) -> Vec<f64> {
    (2..e.len()).map(|i| (e[i] / e[i - 1]).ln() / (e[i - 1] / e[i - 2]).ln()).collect()
}

fn main() {
    // (a) c50 and lambda-hat
    let (r, k) = secant(f, 1.0, 2.0, 1e-12, 50);
    println!("c50        = {:.12} after {} steps", r, k);
    let (r, k) = secant(score, 2.0, 2.5, 1e-12, 50);
    println!("lambda-hat = {:.12} after {} steps", r, k);
    // (b) errors and order estimates on the score
    let lhat = 1.2184695258301268; // reference: Brent (uniroot/brentq) with tol 1e-15
    let mut xs = vec![2.0, 2.5];
    for _ in 0..9 {
        let m = xs.len();
        let (a, b) = (xs[m - 2], xs[m - 1]);
        xs.push(b - score(b) * (b - a) / (score(b) - score(a)));
    }
    let mut xn = vec![2.0];
    for k in 0..5 {
        xn.push(xn[k] - score(xn[k]) / dscore(xn[k]));
    }
    // errors, dropping those at rounding level
    let es: Vec<f64> = xs.iter().map(|x| (x - lhat).abs()).filter(|e| *e > 1e-14).collect();
    let en: Vec<f64> = xn.iter().map(|x| (x - lhat).abs()).filter(|e| *e > 1e-14).collect();
    for e in [&es, &en] {
        println!("{:?}", e.iter().map(|v| format!("{:.2e}", v)).collect::<Vec<_>>());
        println!("{:?}", order_est(e).iter().map(|v| format!("{:.2}", v)).collect::<Vec<_>>());
    }
}
