// C5: the repressilator with ode_solvers: Dopri5 at several tolerances, then fixed-step Rk4
use ode_solvers::{Dopri5, OutputType, Rk4, System, Vector6};

type State = Vector6<f64>;
const ALPHA: f64 = 216.0; const ALPHA0: f64 = 0.216; const BETA: f64 = 5.0; const N: f64 = 2.0;

struct Repressilator;
impl System<f64, State> for Repressilator {
    fn system(&self, _t: f64, y: &State, dy: &mut State) {
        let rep = [y[5], y[3], y[4]]; // m1 repressed by p3, m2 by p1, m3 by p2
        for i in 0..3 {
            dy[i] = -y[i] + ALPHA / (1.0 + rep[i].powf(N)) + ALPHA0;
            dy[i + 3] = -BETA * (y[i + 3] - y[i]);
        }
    }
}

fn main() {
    let y0 = State::from_column_slice(&[1.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
    // (a) Dopri5 with tolerances 1e-3 ... 1e-9
    for k in 3..10 {
        let tol = 10f64.powi(-k);
        let mut s = Dopri5::new(Repressilator, 0.0, 50.0, 1.0, y0, tol, tol);
        s.set_output(OutputType::Sparse); // last output = the step that ends at t = 50
        let st = s.integrate().unwrap();
        println!("tol {:e}: p1(50) = {:.8}, steps {}, f evaluations {}",
                 tol, s.y_out().last().unwrap()[3], st.accepted_steps, st.num_eval);
    }
    // (b) fixed-step RK4: Rk4::new(system, t0, y0, t_end, h)
    for h in [1.0, 0.2, 0.01] {
        let mut s = Rk4::new(Repressilator, 0.0, y0, 50.0, h);
        s.integrate().unwrap();
        println!("rk4, h = {}: p1(50) = {:.8e}", h, s.y_out().last().unwrap()[3]);
    }
}
