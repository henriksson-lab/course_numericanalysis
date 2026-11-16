// Q5: Newton on a double root, f(x) = (x - 1)^2, converges only linearly.
fn main() {
    let f = |x: f64| (x - 1.0).powi(2);
    let df = |x: f64| 2.0 * (x - 1.0);
    let mut x: f64 = 2.0;
    for k in 1..=6 {
        x -= f(x) / df(x);
        println!("{} {:.6} {:.6}", k, x, x - 1.0);
    }
}
