// C2: cost per step versus accuracy (scalar versions of the C1 steppers)
type Step = fn(&dyn Fn(f64) -> f64, f64, f64) -> f64; // autonomous: f depends on N only

fn euler_step(f: &dyn Fn(f64) -> f64, y: f64, h: f64) -> f64 {
    y + h * f(y)
}
fn heun_step(f: &dyn Fn(f64) -> f64, y: f64, h: f64) -> f64 {
    let k1 = f(y);
    let k2 = f(y + h * k1);
    y + h / 2.0 * (k1 + k2)
}
fn rk4_step(f: &dyn Fn(f64) -> f64, y: f64, h: f64) -> f64 {
    let k1 = f(y);
    let k2 = f(y + h / 2.0 * k1);
    let k3 = f(y + h / 2.0 * k2);
    let k4 = f(y + h * k3);
    y + h / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
}

fn main() {
    let f = |n: f64| 0.8 * n * (1.0 - n / 1.2);
    let exact = |t: f64| 1.2 / (1.0 + (1.2 / 0.05 - 1.0) * (-0.8 * t).exp());
    let methods: [(&str, Step, usize); 3] =
        [("Euler", euler_step, 1), ("Heun", heun_step, 2), ("RK4", rk4_step, 4)];
    for (name, step, stages) in methods {
        for k in 0..14 {
            let n = 10 * 2usize.pow(k); // h = 1, 1/2, ..., 1/8192
            let h = 10.0 / n as f64;
            let mut y = 0.05;
            for _ in 0..n {
                y = step(&f, y, h);
            }
            let e = (y - exact(10.0)).abs();
            if e < 1e-6 {
                println!("{}: {} steps, {} evaluations, error {:.3e}", name, n, n * stages, e);
                break;
            }
        }
    }
}
