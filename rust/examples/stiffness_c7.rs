// C7: backward Euler damps a 24-hour oscillation; an explicit library solver does not.
use nalgebra::{DMatrix, DVector};
use ode_solvers::{Dopri5, System, Vector2};
use std::f64::consts::PI;

const W: f64 = 2.0 * PI / 24.0;                       // one cycle per 24 h

fn amp(x: f64, v: f64) -> f64 { (x * x + (v / W).powi(2)).sqrt() } // exact: 1 for ever

struct Osc;
impl System<f64, Vector2<f64>> for Osc {
    fn system(&self, _t: f64, y: &Vector2<f64>, dy: &mut Vector2<f64>) {
        dy[0] = y[1]; dy[1] = -W * W * y[0];
    }
}

fn main() {
    let a = DMatrix::from_row_slice(2, 2, &[0.0, 1.0, -W * W, 0.0]);
    let id = DMatrix::<f64>::identity(2, 2);
    for h in [1.0, 0.5, 0.1, 0.01] {
        // linear ODE: one Newton iteration is exact, so solve (I - hA) y_{n+1} = y_n directly
        let lu = (&id - h * &a).lu();
        let mut y = DVector::from_vec(vec![1.0, 0.0]);
        for _ in 0..(240.0 / h).round() as usize { y = lu.solve(&y).unwrap(); }
        println!("backward Euler h = {:4.2}: amplitude after 10 days {:.6}", h, amp(y[0], y[1]));
    }
    let mut s = Dopri5::new(Osc, 0.0, 240.0, 240.0, Vector2::new(1.0, 0.0), 1e-8, 1e-10);
    let st = s.integrate().unwrap();
    let y = s.y_out().last().unwrap();
    println!("Dopri5: amplitude {:.6}, steps {}", amp(y[0], y[1]), st.accepted_steps);
}
