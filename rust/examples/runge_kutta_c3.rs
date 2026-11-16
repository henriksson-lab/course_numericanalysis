// C3: predator-prey; drift of the conserved quantity V for fixed steps and for Dopri5
use ode_solvers::{Dopri5, OutputType, System, Vector2};

const A: f64 = 1.0; const B: f64 = 0.1; const C: f64 = 0.075; const D: f64 = 1.5; // per year
type State = Vector2<f64>;

fn lv(y: &State) -> State {
    State::new(A * y[0] - B * y[0] * y[1], C * y[0] * y[1] - D * y[1])
}
fn v(y: &State) -> f64 {
    C * y[0] - D * y[0].ln() + B * y[1] - A * y[1].ln()
}
fn drift(ys: &[State]) -> f64 {
    ys.iter().map(|y| (v(y) - v(&ys[0])).abs()).fold(0.0, f64::max)
}

// nalgebra vectors support y + k * h directly
fn euler_step(y: &State, h: f64) -> State { y + lv(y) * h }
fn heun_step(y: &State, h: f64) -> State {
    let k1 = lv(y);
    let k2 = lv(&(y + k1 * h));
    y + (k1 + k2) * (h / 2.0)
}
fn rk4_step(y: &State, h: f64) -> State {
    let k1 = lv(y);
    let k2 = lv(&(y + k1 * (h / 2.0)));
    let k3 = lv(&(y + k2 * (h / 2.0)));
    let k4 = lv(&(y + k3 * h));
    y + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (h / 6.0)
}

struct LotkaVolterra;
impl System<f64, State> for LotkaVolterra {
    fn system(&self, _t: f64, y: &State, dy: &mut State) { *dy = lv(y); }
}

fn main() {
    let methods: [(&str, fn(&State, f64) -> State, f64); 3] =
        [("Euler", euler_step, 0.01), ("Heun", heun_step, 0.01), ("RK4", rk4_step, 0.1)];
    for (name, step, h) in methods {
        let mut ys = vec![State::new(10.0, 5.0)];
        for _ in 0..(50.0 / h).round() as usize {
            ys.push(step(ys.last().unwrap(), h));
        }
        println!("{}: drift {:.4e}", name, drift(&ys));
    }
    for tol in [1e-3, 1e-6] {
        let mut s = Dopri5::new(LotkaVolterra, 0.0, 50.0, 1.0, State::new(10.0, 5.0), tol, tol);
        s.set_output(OutputType::Sparse); // keep the solver's own steps
        let stats = s.integrate().unwrap();
        println!("Dopri5, tol = {:e}: drift {:.3e}, steps {}", tol, drift(s.y_out()), stats.accepted_steps);
    }
}
