// Background: one hand-written RK4 step, then the library solver (ode_solvers::Dopri5)
use ode_solvers::continuous_output_model::ContinuousOutputModel;
use ode_solvers::{Dopri5, System, Vector1};

fn f(_t: f64, n: f64) -> f64 {
    0.8 * n * (1.0 - n / 1.2) // logistic, per hour
}

fn rk4_step(f: fn(f64, f64) -> f64, t: f64, y: f64, h: f64) -> f64 {
    let k1 = f(t, y);
    let k2 = f(t + h / 2.0, y + h / 2.0 * k1);
    let k3 = f(t + h / 2.0, y + h / 2.0 * k2);
    let k4 = f(t + h, y + h * k3);
    y + h / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
}

// ode_solvers wants the ODE as a type implementing the System trait
struct Logistic;
impl System<f64, Vector1<f64>> for Logistic {
    fn system(&self, t: f64, y: &Vector1<f64>, dy: &mut Vector1<f64>) {
        dy[0] = f(t, y[0]);
    }
}

fn main() {
    println!("{}", rk4_step(f, 0.0, 0.05, 1.0)); // one RK4 step of h = 1

    // Dopri5::new(system, t0, t_end, output spacing, y0, rtol, atol)
    let mut solver = Dopri5::new(Logistic, 0.0, 10.0, 1.0, Vector1::new(0.05), 1e-8, 1e-10);
    let mut dense = ContinuousOutputModel::default(); // the solution at any t
    solver.integrate_with_continuous_output_model(&mut dense).unwrap();
    for t in [0.0, 1.0, 10.0] {
        println!("t = {:>2}: N = {:.7}", t, dense.evaluate(t).unwrap()[0]);
    }
}
