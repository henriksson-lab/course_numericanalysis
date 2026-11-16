// rust/examples/phase_plane_c3.rs  —  cargo run --example phase_plane_c3
// C3: where does the cell end up? Simulate from a 10 x 10 grid of starts, for the
// symmetric switch and for the biased one (alpha_v = 2.9). Steady states of the
// biased switch by a hand-written Newton (no root-finding crate in our list).
// ode_solvers uses its own nalgebra version, so Newton uses nalgebra's types
// (prefix na::) and we convert with Vector2::new(x[0], x[1]).
use nalgebra as na;
use ode_solvers::{Dopri5, System, Vector2};

struct Toggle { au: f64, av: f64 }
impl System<f64, Vector2<f64>> for Toggle {
    fn system(&self, _t: f64, x: &Vector2<f64>, dx: &mut Vector2<f64>) {
        dx[0] = self.au / (1.0 + x[1] * x[1]) - x[0];
        dx[1] = self.av / (1.0 + x[0] * x[0]) - x[1];
    }
}

fn newton(au: f64, av: f64, u0: f64, v0: f64) -> Vector2<f64> {
    let mut x = na::Vector2::new(u0, v0);
    for _ in 0..50 {
        let (u, v) = (x[0], x[1]);
        let fx = na::Vector2::new(au / (1.0 + v * v) - u, av / (1.0 + u * u) - v);
        let j = na::Matrix2::new(-1.0, -2.0 * au * v / (1.0 + v * v).powi(2),
                                 -2.0 * av * u / (1.0 + u * u).powi(2), -1.0);
        x += j.lu().solve(&(-fx)).expect("singular Jacobian");
    }
    Vector2::new(x[0], x[1])
}

fn which_state(x: &Vector2<f64>, states: &[(&str, Vector2<f64>)], tol: f64) -> String {
    let (name, d) = states.iter().map(|(k, s)| (*k, (x - s).norm()))
        .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap()).unwrap();
    if d < tol { name.to_string() } else { "not settled".to_string() }
}

fn fates(au: f64, av: f64, states: &[(&str, Vector2<f64>)]) {
    let mut counts = vec![0; states.len() + 1];
    for i in 0..10 {
        for k in 0..10 {
            let x0 = Vector2::new(0.2 + 0.4 * i as f64, 0.2 + 0.4 * k as f64);
            let mut s = Dopri5::new(Toggle { au, av }, 0.0, 100.0, 1.0, x0, 1e-10, 1e-12);
            s.integrate().expect("integration failed");
            let name = which_state(s.y_out().last().unwrap(), states, 1e-3);
            let idx = states.iter().position(|(k, _)| *k == name).unwrap_or(states.len());
            counts[idx] += 1;
        }
    }
    let names: Vec<&str> = states.iter().map(|(k, _)| *k).chain(["not settled"]).collect();
    println!("{:?}", names.iter().zip(&counts).collect::<Vec<_>>());
}

fn main() {
    let sym = [("U on", Vector2::new(2.618034, 0.381966)),
               ("middle", Vector2::new(1.213412, 1.213412)),
               ("V on", Vector2::new(0.381966, 2.618034))];
    fates(3.0, 3.0, &sym);
    let (au, av) = (3.0, 2.9); // biased switch
    let biased = [("U on", newton(au, av, 3.0, 0.3)),
                  ("middle", newton(au, av, 1.1, 1.3)),
                  ("V on", newton(au, av, 0.4, 2.5))];
    for (k, s) in &biased { println!("{:<7} {:.7} {:.7}", k, s[0], s[1]); }
    fates(au, av, &biased);
}
