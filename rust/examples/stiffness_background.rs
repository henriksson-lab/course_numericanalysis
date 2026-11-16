// Background: the enzyme model E + S <-> C -> E + P with ode_solvers::Dopri5 (explicit RK45).
// ode_solvers has no stiff (implicit) solver: an explicit method works here, but needs many steps.
use ode_solvers::continuous_output_model::ContinuousOutputModel;
use ode_solvers::{Dopri5, System, Vector4};

const KF: f64 = 100.0; const KR: f64 = 100.0; const KCAT: f64 = 1.0; // per uM per s, per s, per s

struct Enzyme;
impl System<f64, Vector4<f64>> for Enzyme {
    fn system(&self, _t: f64, y: &Vector4<f64>, dy: &mut Vector4<f64>) {
        let (s, e, c) = (y[0], y[1], y[2]);
        let bind = KF * e * s - KR * c; // E + S <-> C  (net binding)
        let turn = KCAT * c;            // C -> E + P   (catalysis)
        dy[0] = -bind; dy[1] = -bind + turn; dy[2] = bind - turn; dy[3] = turn;
    }
}

fn main() {
    let y0 = Vector4::new(10.0, 1.0, 0.0, 0.0); // uM: S, E, C, P
    // Dopri5::new(system, t0, t_end, output spacing, y0, rtol, atol)
    let mut s = Dopri5::new(Enzyme, 0.0, 30.0, 30.0, y0, 1e-8, 1e-10);
    let mut dense = ContinuousOutputModel::default(); // the solution at any t
    let st = s.integrate_with_continuous_output_model(&mut dense).unwrap();
    for t in [10.0, 30.0] {
        let y = dense.evaluate(t).unwrap();
        println!("t = {:>2}  S = {:.6e}  E = {:.6}  C = {:.6e}  P = {:.6}", t, y[0], y[1], y[2], y[3]);
    }
    println!("{} steps, {} f evaluations", st.accepted_steps, st.num_eval);
}
