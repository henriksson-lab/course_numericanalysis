// rust/examples/phase_plane_c1.rs  —  cargo run --example phase_plane_c1
// C1: trajectories of the toggle switch from seven starts; where are they at t = 30 and t = 60?
// (No plotting here: write ts/xs to a file and plot with any tool if you like.)
use ode_solvers::{Dopri5, System, Vector2};

const ALPHA: f64 = 3.0;
const N: f64 = 2.0;

struct Toggle;
impl System<f64, Vector2<f64>> for Toggle {
    fn system(&self, _t: f64, x: &Vector2<f64>, dx: &mut Vector2<f64>) {
        let (u, v) = (x[0], x[1]);
        dx[0] = ALPHA / (1.0 + v.powf(N)) - u;
        dx[1] = ALPHA / (1.0 + u.powf(N)) - v;
    }
}

const STATES: [(&str, [f64; 2]); 3] = [("U on", [2.618034, 0.381966]),
                                       ("middle", [1.213412, 1.213412]),
                                       ("V on", [0.381966, 2.618034])];

// name of the nearest state, if closer than tol
fn which_state(x: &Vector2<f64>, states: &[(&str, [f64; 2])], tol: f64) -> String {
    let (name, d) = states.iter()
        .map(|(k, s)| (*k, ((x[0] - s[0]).powi(2) + (x[1] - s[1]).powi(2)).sqrt()))
        .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap()).unwrap();
    if d < tol { name.to_string() } else { "not settled".to_string() }
}

fn main() {
    let starts = [(0.5, 0.3), (0.3, 0.5), (3.3, 3.0), (3.0, 3.3), (0.1, 1.8), (1.8, 0.1), (2.0, 2.0)];
    for (u0, v0) in starts {
        // output every 0.05 time units, as in R and Python
        let mut s = Dopri5::new(Toggle, 0.0, 60.0, 0.05, Vector2::new(u0, v0), 1e-10, 1e-12);
        s.integrate().expect("integration failed");
        let xs = s.y_out();
        let (x30, x60) = (&xs[600], xs.last().unwrap()); // t = 30 is entry 600
        println!("({:.1}, {:.1})  t = 30: {:<11}  t = 60: ({:.6}, {:.6}) {}", u0, v0,
                 which_state(x30, &STATES, 1e-3), x60[0], x60[1], which_state(x60, &STATES, 1e-3));
    }
}
