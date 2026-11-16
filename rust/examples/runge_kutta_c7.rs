// C7: events by hand: ode_solvers has no event location, so we bracket the sign
// change of g on the dense output and refine it by bisection
use ode_solvers::continuous_output_model::ContinuousOutputModel;
use ode_solvers::{Dopri5, System, Vector2, Vector3};

struct Sir;
impl System<f64, Vector3<f64>> for Sir {
    fn system(&self, _t: f64, y: &Vector3<f64>, dy: &mut Vector3<f64>) {
        let (beta, gamma, n) = (0.3, 0.1, 10000.0);
        let (s, i) = (y[0], y[1]);
        *dy = Vector3::new(-beta * s * i / n, beta * s * i / n - gamma * i, gamma * i);
    }
}
struct LotkaVolterra;
impl System<f64, Vector2<f64>> for LotkaVolterra {
    fn system(&self, _t: f64, y: &Vector2<f64>, dy: &mut Vector2<f64>) {
        let (a, b, c, d) = (1.0, 0.1, 0.075, 1.5);
        *dy = Vector2::new(a * y[0] - b * y[0] * y[1], c * y[0] * y[1] - d * y[1]);
    }
}

// first root of g(t) on [t0, t1]: scan with spacing dt, then bisect
fn first_root(g: &dyn Fn(f64) -> f64, t0: f64, t1: f64, dt: f64) -> Option<f64> {
    let mut lo = t0;
    while lo < t1 {
        let hi = (lo + dt).min(t1);
        if g(lo).signum() != g(hi).signum() {
            let (mut a, mut b) = (lo, hi);
            for _ in 0..60 {
                let m = 0.5 * (a + b);
                if g(a).signum() == g(m).signum() { a = m } else { b = m }
            }
            return Some(0.5 * (a + b));
        }
        lo = hi;
    }
    None
}

fn main() {
    // (a) SIR: the peak as an event, g = beta S / N - gamma
    let mut dense = ContinuousOutputModel::default();
    let mut s = Dopri5::new(Sir, 0.0, 160.0, 1.0, Vector3::new(9990.0, 10.0, 0.0), 1e-10, 1e-8);
    s.integrate_with_continuous_output_model(&mut dense).unwrap();
    let g = |t: f64| 0.3 * dense.evaluate(t).unwrap()[0] / 10000.0 - 0.1;
    let tp = first_root(&g, 0.0, 160.0, 0.1).unwrap();
    let y = dense.evaluate(tp).unwrap();
    println!("peak at t = {:.8}: S = {:.6}, I = {:.6}, R = {:.6}", tp, y[0], y[1], y[2]);

    // (b) predator-prey: first time the lynx exceed 20
    let mut dense = ContinuousOutputModel::default();
    let mut s = Dopri5::new(LotkaVolterra, 0.0, 15.0, 1.0, Vector2::new(10.0, 5.0), 1e-10, 1e-10);
    s.integrate_with_continuous_output_model(&mut dense).unwrap();
    let g = |t: f64| dense.evaluate(t).unwrap()[1] - 20.0;
    let t20 = first_root(&g, 0.0, 15.0, 0.1).unwrap();
    println!("lynx reach 20 at t = {:.8}, prey x = {:.6}", t20, dense.evaluate(t20).unwrap()[0]);

    // (c) only looking at t = 0, 2.5, ..., 15
    let ys: Vec<f64> = (0..7).map(|k| dense.evaluate(2.5 * k as f64).unwrap()[1]).collect();
    println!("{:.2?}", ys);
    println!("found with spacing 2.5: {:?}", first_root(&g, 0.0, 15.0, 2.5));
}
