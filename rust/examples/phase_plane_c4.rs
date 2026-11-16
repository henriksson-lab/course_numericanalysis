// rust/examples/phase_plane_c4.rs  —  cargo run --example phase_plane_c4
// C4: the separatrix of the biased switch (alpha_v = 2.9), from two runs backwards in
// time out of the saddle along its stable eigenvector; then check it with forward runs.
// Backwards in time = forwards with the sign of the right-hand side flipped.
use nalgebra as na;
use ode_solvers::{Dopri5, System, Vector2};

const AU: f64 = 3.0;
const AV: f64 = 2.9;

fn rhs(x: &Vector2<f64>) -> Vector2<f64> {
    Vector2::new(AU / (1.0 + x[1] * x[1]) - x[0], AV / (1.0 + x[0] * x[0]) - x[1])
}

struct Toggle { dir: f64 } // dir = 1: forwards in time, dir = -1: backwards
impl System<f64, Vector2<f64>> for Toggle {
    fn system(&self, _t: f64, x: &Vector2<f64>, dx: &mut Vector2<f64>) {
        *dx = rhs(x) * self.dir;
    }
}

fn run(x0: Vector2<f64>, t_end: f64, dir: f64) -> Vec<Vector2<f64>> {
    let mut s = Dopri5::new(Toggle { dir }, 0.0, t_end, 0.01, x0, 1e-10, 1e-12);
    s.integrate().expect("integration failed");
    s.y_out().clone()
}

fn newton(mut x: Vector2<f64>) -> Vector2<f64> {
    for _ in 0..50 {
        let j = fd_jacobian(&x);
        let fx = rhs(&x);
        let d = j.lu().solve(&na::Vector2::new(-fx[0], -fx[1])).expect("singular Jacobian");
        x += Vector2::new(d[0], d[1]);
    }
    x
}

fn fd_jacobian(x: &Vector2<f64>) -> na::Matrix2<f64> {
    let h = f64::EPSILON.sqrt();
    let fx = rhs(x);
    let mut jm = na::Matrix2::zeros();
    for j in 0..2 {
        let hj = h * x[j].abs().max(1.0);
        let mut xh = *x;
        xh[j] += hj;
        let col = (rhs(&xh) - fx) / hj;
        jm[(0, j)] = col[0];
        jm[(1, j)] = col[1];
    }
    jm
}

fn which_state(x: &Vector2<f64>, states: &[(&str, Vector2<f64>)]) -> String {
    let (name, d) = states.iter().map(|(k, s)| (*k, (x - s).norm()))
        .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap()).unwrap();
    if d < 1e-3 { name.to_string() } else { "not settled".to_string() }
}

fn main() {
    let states = [("U on", newton(Vector2::new(3.0, 0.3))),
                  ("middle", newton(Vector2::new(1.1, 1.3))),
                  ("V on", newton(Vector2::new(0.4, 2.5)))];
    let xs = states[1].1;
    let j = fd_jacobian(&xs);
    // the negative eigenvalue; real here (det J < 0), so take the real part
    let lam = j.complex_eigenvalues().iter().map(|z| z.re).fold(f64::INFINITY, f64::min);
    // eigenvector of a 2x2 by hand: (J - lam I) v = 0 is solved by v = (J12, lam - J11)
    let vs = Vector2::new(j[(0, 1)], lam - j[(0, 0)]).normalize();
    let mut sep: Vec<Vector2<f64>> = Vec::new();
    for sgn in [1.0, -1.0] {
        sep.extend(run(xs + vs * (sgn * 1e-6), 7.5, -1.0)); // 7.5 time units backwards
    }
    sep.retain(|p| p[0] >= 0.0 && p[0] <= 4.0);
    sep.sort_by(|p, q| p[0].partial_cmp(&q[0]).unwrap()); // sort by u for interpolation
    for u0 in [0.1, 0.5, 1.0, 2.0, 3.0] {
        let k = sep.iter().position(|p| p[0] >= u0).unwrap(); // linear interpolation
        let (p, q) = (sep[k - 1], sep[k]);
        let vsep = p[1] + (u0 - p[0]) / (q[0] - p[0]) * (q[1] - p[1]);
        let below = *run(Vector2::new(u0, vsep - 0.01), 200.0, 1.0).last().unwrap();
        let above = *run(Vector2::new(u0, vsep + 0.01), 200.0, 1.0).last().unwrap();
        println!("u = {:.1}  v_sep = {:.4}  below -> {:<5} above -> {}", u0, vsep,
                 which_state(&below, &states), which_state(&above, &states));
    }
}
