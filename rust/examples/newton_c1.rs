// C1: check the derivative with a forward difference, then watch Newton's iterates.
fn f(c: f64) -> f64 {
    0.3 * c / (0.1 + c) + 0.7 * c.powi(2) / (25.0 + c.powi(2)) - 0.5
}
fn df(c: f64) -> f64 {
    0.03 / (0.1 + c).powi(2) + 35.0 * c / (25.0 + c.powi(2)).powi(2)
}

fn main() {
    // (a) derivative check with a forward difference
    let h = 1e-6;
    for c in [1.0, 3.0, 10.0] {
        let fd = (f(c + h) - f(c)) / h;
        println!("{:4.0}  {:.10}  {:.10}  {:.1e}", c, df(c), fd, (fd - df(c)).abs() / df(c));
    }
    // (b) Newton, printing every iterate and the step size
    for x0 in [1.0, 5.0] {
        println!("start {}", x0);
        let mut x: f64 = x0;
        for k in 1..=10 {
            let step = f(x) / df(x);
            x -= step;
            println!("{:2}  {:.15}  {:.1e}", k, x, step.abs());
            if step.abs() <= 1e-12 * x.abs() {
                break;
            }
        }
    }
}
