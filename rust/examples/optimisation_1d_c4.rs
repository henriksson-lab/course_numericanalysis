// C4: Newton for minimisation of the decay-rate SSE, with derivative checks and
// the two safety checks (f'' > 0, f decreases).
const T: [f64; 7] = [0.0, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0];
const Y: [f64; 7] = [1.000, 0.804, 0.684, 0.484, 0.217, 0.094, 0.083];

fn sse(k: f64) -> f64 {
    T.iter().zip(Y.iter()).map(|(t, y)| (y - (-k * t).exp()).powi(2)).sum()
}
fn d1(k: f64) -> f64 { // SSE'
    T.iter().zip(Y.iter()).map(|(t, y)| 2.0 * (y - (-k * t).exp()) * t * (-k * t).exp()).sum()
}
fn d2(k: f64) -> f64 { // SSE''
    T.iter().zip(Y.iter())
        .map(|(t, y)| 2.0 * t * t * (-k * t).exp() * (2.0 * (-k * t).exp() - y)).sum()
}

// Some((xmin, iterations)) on success, None (with a message) when a check fails
fn newton_min(f: fn(f64) -> f64, d1: fn(f64) -> f64, d2: fn(f64) -> f64,
              x0: f64, tol: f64, maxit: usize) -> Option<(f64, usize)> {
    let mut x = x0;
    for k in 1..=maxit {
        if d2(x) <= 0.0 { // local parabola opens downwards
            println!("stopped: f'' <= 0 at x = {:.4}", x);
            return None;
        }
        let xnew = x - d1(x) / d2(x);
        if (xnew - x).abs() <= tol * xnew.abs() {
            return Some((xnew, k));
        }
        if !f(xnew).is_finite() || f(xnew) > f(x) { // the step went uphill
            println!("stopped: f increases, step from {:.4} to {:.4}", x, xnew);
            return None;
        }
        x = xnew;
    }
    println!("no convergence in {} steps", maxit);
    None
}

fn main() {
    let h = 1e-5; // (a) check against central differences
    println!("{:.6} {:.6}", d1(0.3), (sse(0.3 + h) - sse(0.3 - h)) / (2.0 * h));
    println!("{:.5} {:.5}", d2(0.3), (d1(0.3 + h) - d1(0.3 - h)) / (2.0 * h));
    if let Some((xmin, it)) = newton_min(sse, d1, d2, 0.2, 1e-12, 50) {
        println!("{:.10} {}", xmin, it);
    }
    println!("{:?}", newton_min(sse, d1, d2, 0.7, 1e-12, 50));
    println!("{:?}", newton_min(sse, d1, d2, 1.0, 1e-12, 50));
}
